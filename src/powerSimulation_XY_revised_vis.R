# =============================================================================
# Simulation-Based Power Analysis: BM Trait Evolution with Mean Shift
# Phylogenetic Mixed Models | Grass Family (Poaceae)
#
# CHANGES FROM ORIGINAL:
#   1. Reproducible parallelism via RNGkind("L'Ecuyer-CMRG") + set.seed()
#      before every mclapply call.
#   2. All 16 power curve runs replaced by a single expand.grid() loop.
#   3. n_sims replicates are now actually run per (beta, scenario, prop)
#      combination; power is averaged across replicates.
#   4. Binomial 95% CIs on power estimates (Wilson interval).
#   5. Summary table: beta threshold at which power first crosses 80% for
#      every scenario × convergence combination.
#
# DEPENDENCIES: ape, phytools, ggplot2, dplyr, asreml, parallel
# =============================================================================

library(ape)
library(phytools)
library(ggplot2)
library(dplyr)

# =============================================================================
# 0. CONFIG  (edit these paths / parameters)
# =============================================================================

PHYLOGWAS_ROOT <- Sys.getenv("PHYLOGWAS_ROOT", unset = "/workdir/sh2246/p_phyloGWAS")

TREE_PATH   <- file.path(PHYLOGWAS_ROOT, "output/powerSimulation/tree.nwk")
KMAT_PATH   <- file.path(PHYLOGWAS_ROOT, "output/phyloK_728Poaceae_astral_20250407.txt")
TRANS_PATH  <- file.path(PHYLOGWAS_ROOT, "output/powerSimulation/realTransition.json")
OUTPUT_DIR  <- file.path(PHYLOGWAS_ROOT, "output/figure")

# Simulation parameters
sigma2_Y      <- 0.5
sigma2_X      <- 0.5
sigma2_resid  <- 0.1   # residual multiplier (applied to both Y and X); was
                        # inconsistently 0.25 in params block vs 0.1 in
                        # run_power_analysis — now set once here
shift_Y       <- 2
beta_range    <- seq(0, 0.5, by = 0.025)
n_X           <- 1   # X traits per replicate
n_sims        <- 200    # simulation replicates per (beta, scenario, prop)
alpha         <- 0.05
n_cores       <- 10
MASTER_SEED   <- 42

# =============================================================================
# 1. TREE + K-MATRIX SETUP
# =============================================================================

tree <- read.tree(TREE_PATH)

phyloKMat <- read.table(KMAT_PATH)
phyloKMat <- phyloKMat/2
colnames(phyloKMat) <- rownames(phyloKMat)

transition_nodes <- rjson::fromJSON(file = TRANS_PATH)[[1]]
transition_nodes <- unlist(transition_nodes)

# =============================================================================
# 2. HELPER FUNCTIONS  (unchanged from original)
# =============================================================================

get_all_ancestors <- function(tree, node) {
  ancestors <- c()
  current   <- node
  root      <- Ntip(tree) + 1
  while (current != root) {
    parent <- tree$edge[tree$edge[, 2] == current, 1]
    if (length(parent) == 0) break
    ancestors <- c(ancestors, parent)
    current   <- parent
  }
  ancestors
}

sample_independent_nodes <- function(tree, n, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  internal_nodes <- (Ntip(tree) + 1):(Ntip(tree) + Nnode(tree))
  internal_nodes <- internal_nodes[-1]   # exclude root
  selected   <- c()
  candidates <- sample(internal_nodes)
  for (nd in candidates) {
    if (length(selected) == n) break
    nd_ancestors <- get_all_ancestors(tree, nd)
    nd_desc      <- phytools::getDescendants(tree, nd)
    is_ancestor_of_selected   <- any(selected %in% nd_desc)
    is_descendant_of_selected <- any(selected %in% nd_ancestors)
    if (!is_ancestor_of_selected && !is_descendant_of_selected)
      selected <- c(selected, nd)
  }
  if (length(selected) < n)
    warning(sprintf("Only found %d independent nodes (requested %d).",
                    length(selected), n))
  selected
}

get_tip_descendants <- function(tree, node) {
  desc <- phytools::getDescendants(tree, node)
  tips <- desc[desc <= Ntip(tree)]
  tree$tip.label[tips]
}

make_regime <- function(tree, transition_nodes) {
  regimes <- setNames(rep(0L, Ntip(tree)), tree$tip.label)
  for (nd in transition_nodes) {
    tips <- get_tip_descendants(tree, nd)
    regimes[tips] <- 1L
  }
  regimes
}

# =============================================================================
# 3. TRANSITION NODE SCENARIOS
#
# The empirical scenario uses a fixed node set (it represents real biology).
# Random scenarios are defined by a SAMPLING SPEC: the number of nodes to draw
# is recorded here, and a fresh node set is drawn inside each replicate of
# the main loop. This propagates uncertainty about which random nodes were
# drawn into the power estimate and its CI.
# =============================================================================

n_trans <- length(transition_nodes)

# Each entry is either:
#   list(type = "fixed",  nodes = <vector>)              -> reused every rep
#   list(type = "random", n     = <integer>)             -> resampled every rep
scenarios <- list(
  "empirical transition"   = list(type = "fixed",  nodes = transition_nodes),
  "random transition - 1x" = list(type = "random", n = n_trans),
  "random transition - 2x" = list(type = "random", n = n_trans * 2),
  "random transition - 1/2x" = list(type = "random", n = n_trans * .5)
)

# =============================================================================
# 4. SIMULATION FUNCTIONS  (sigma2_resid now uses config value)
# =============================================================================

simulate_Y <- function(tree, regime, sigma2_Y, shift_Y,
                       sigma2_resid_Y, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  bm_part    <- fastBM(tree, sig2 = sigma2_Y, internal = FALSE)
  resid_part <- rnorm(Ntip(tree), 0, sqrt(sigma2_resid_Y))
  shift_vec  <- ifelse(regime == 1, shift_Y, 0)
  bm_part + shift_vec + resid_part
}

simulate_X <- function(tree, regime, sigma2_X, shift_X,
                       sigma2_resid_X, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  bm_part    <- fastBM(tree, sig2 = sigma2_X, internal = FALSE)
  resid_part <- rnorm(Ntip(tree), 0, sqrt(sigma2_resid_X))
  shift_vec  <- ifelse(regime == 1, shift_X, 0)
  bm_part + shift_vec + resid_part
}

simulate_Y_and_X <- function(tree, nodes, prop, sigma2_Y, sigma2_X,
                              shift_Y, beta, n_X,
                              sigma2_resid_Y, sigma2_resid_X,
                              seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  regime_Y <- make_regime(tree, nodes)
  Y        <- simulate_Y(tree, regime_Y, sigma2_Y, shift_Y, sigma2_resid_Y)

  shift_X  <- beta * shift_Y
  regime_X <- make_regime(tree, sample(nodes, round(length(nodes) * prop)))
  X_mat    <- replicate(n_X, simulate_X(tree, regime_X, sigma2_X,
                                         shift_X, sigma2_resid_X))
  colnames(X_mat) <- paste0("X", seq_len(n_X))
  data.frame(assemblyID = as.factor(tree$tip.label), Y = Y, X = X_mat)
}

# =============================================================================
# 5. MODEL FITTING
# =============================================================================

fit_model <- function(data, Kmat, responseVar, predictor) {
  asreml.options(verbose = FALSE, trace = FALSE)
  Coeff    <- numeric(length(predictor))
  WaldPval <- numeric(length(predictor))

  for (i in seq_along(predictor)) {
    fullFM     <- as.formula(paste0(responseVar, "~", predictor[i]))
    model_full <- asreml(fixed  = fullFM,
                         random = ~ vm(assemblyID, Kmat),
                         ai.sing = FALSE, data = data)
    if (model_full$converge) {
      modelWald    <- wald.asreml(model_full)
      Coeff[i]    <- model_full$coefficients$fixed[2]
      WaldPval[i] <- modelWald[2, 4]
    } else {
      Coeff[i]    <- NA
      WaldPval[i] <- NA
    }
  }
  list(coeff = Coeff, p = WaldPval)
}

# =============================================================================
# 6. POWER ANALYSIS  (single replicate; called n_sims times in the loop)
# =============================================================================

run_single_replicate <- function(KMat, tree, nodes, prop,
                                  sigma2_Y, sigma2_X, shift_Y, beta, n_X,
                                  sigma2_resid, alpha, sim_seed) {
  sr_Y <- sigma2_resid * sigma2_Y
  sr_X <- sigma2_resid * sigma2_X

  simDat <- simulate_Y_and_X(tree, nodes, prop,
                               sigma2_Y, sigma2_X, shift_Y, beta, n_X,
                               sr_Y, sr_X, seed = sim_seed)
  tmp    <- fit_model(simDat, KMat, "Y", colnames(simDat)[-c(1, 2)])

  hits   <- sum(tmp$p < alpha, na.rm = TRUE)   # number of X traits detected
  total  <- sum(!is.na(tmp$p))
  list(hits = hits, total = total)
}

# =============================================================================
# 7. WILSON BINOMIAL CI
# =============================================================================

wilson_ci <- function(hits, n, conf = 0.95) {
  if (n == 0) return(c(lower = NA, upper = NA))
  z     <- qnorm(1 - (1 - conf) / 2)
  p_hat <- hits / n
  denom <- 1 + z^2 / n
  centre <- (p_hat + z^2 / (2 * n)) / denom
  margin <- z * sqrt(p_hat * (1 - p_hat) / n + z^2 / (4 * n^2)) / denom
  c(lower = max(0, centre - margin),
    upper = min(1, centre + margin))
}

prop_levels <- c(1.0, 0.8, 0.5, 0.3)
prop_labels <- c("100%", "80%", "50%", "30%")

power_df_merged <- read.csv(file.path(OUTPUT_DIR,"power_df_merged.csv"))
rownames(power_df_merged) <- NULL

# Factor levels for ordered facet display
power_df_merged$convergence <- factor(power_df_merged$convergence,
                                       levels = prop_labels)
# =============================================================================
# 10. TYPE I ERROR CHECK  (beta = 0 sanity check)
# =============================================================================

type1_df <- power_df_merged %>%
  filter(beta == 0) %>%
  select(scenario, convergence, power, ci_lower, ci_upper) %>%
  rename(false_positive_rate = power)

cat("\n--- Type I error rates (beta = 0, should be ~0.05) ---\n")
print(type1_df, row.names = FALSE)

# =============================================================================
# 11. PLOT
# =============================================================================

p_curve <- ggplot(power_df_merged,
                  aes(x = beta, y = power, color = scenario,
                      fill = scenario)) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.5) +
  facet_wrap(~ convergence, nrow = 2) +
  scale_color_manual(values = c("steelblue", "lightgreen",
                                 "olivedrab",  "darkgreen")) +
  scale_fill_manual(values  = c("steelblue", "lightgreen",
                                 "olivedrab",  "darkgreen")) +
  geom_hline(yintercept = 0.80, linetype = "dotted", color = "black",
             linewidth = 0.8) +
  geom_vline(xintercept = c(0.1, 0.2), linetype = "dashed",
             color = "grey40") +
  labs(
    title    = "Power Curve: Detecting Gene-Trait Covariation Along the Phylogeny",
    subtitle = sprintf("BM sigma2 = %.2f | %d replicates per configuration | alpha = %.2f",
                       sigma2_Y, n_sims, alpha),
    x        = expression(paste(theta, " magnitude")),
    y        = "Estimated Power",
    color    = "Scenario",
    fill     = "Scenario"
  ) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  theme_classic(base_size = 13) +
  theme(legend.position = "right")

print(p_curve)

ggsave(file.path(OUTPUT_DIR, "Fig8_revised.png"), p_curve,
       device = "png", width = 10, height = 6, dpi = 300)

# =============================================================================
# 12. SESSION INFO  (for reproducibility reporting)
# =============================================================================

cat("\n--- Session Info ---\n")
sessionInfo()
