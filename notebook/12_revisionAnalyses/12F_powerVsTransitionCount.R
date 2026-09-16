# Follow-up power analysis (not one of the original 5 revision analyses): at a
# FIXED effect size (beta = 0.3) and 100% X/Y regime overlap (the "100%"
# convergence/prop level from src/powerSimulation_XY_revised.R), how many
# independent random transition events are needed on the phylogeny before
# power to detect the effect reaches 0.80?
#
# Sweeps n_trans = 3, 6, 9, ..., 45 (15 values) for the "random transition"
# scenario only -- no empirical/fixed node set needed, so this script doesn't
# touch realTransition.json at all, only the tree + K-matrix. Reuses the exact
# same simulation/fitting logic as src/powerSimulation_XY_revised.R (copied
# in, not sourced, since that script is a full end-to-end script rather than
# a library) so the results are directly comparable to the existing power
# curves.
#
# Usage: N_CORES=<n> Rscript notebook/12_revisionAnalyses/12F_powerVsTransitionCount.R

library(ape)
library(phytools)
library(dplyr)
library(asreml)
library(parallel)

PHYLOGWAS_ROOT <- Sys.getenv("PHYLOGWAS_ROOT", unset = "/workdir/sh2246/p_phyloGWAS")
TREE_PATH  <- Sys.getenv("POWERSIM_TREE_PATH", unset = file.path(PHYLOGWAS_ROOT, "output/powerSimulation/tree.nwk"))
KMAT_PATH  <- Sys.getenv("POWERSIM_KMAT_PATH", unset = file.path(PHYLOGWAS_ROOT, "output/phyloK_728Poaceae_astral_20250407.txt"))
OUTPUT_DIR <- Sys.getenv("POWERSIM_OUTPUT_DIR", unset = file.path(PHYLOGWAS_ROOT, "output/revisionAnalyses/powerSimulation"))
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# --- fixed analysis parameters (this script's whole point: hold these fixed) ---
BETA_FIXED  <- 0.3
PROP_FIXED  <- 1.0     # "100% convergence" -- full overlap between Y's and X's regime nodes
N_TRANS_SEQ <- seq(3, 45, by = 3)   # 3,6,...,45 -- 15 values
POWER_TARGET <- 0.80

sigma2_Y      <- 0.5
sigma2_X      <- 0.5
sigma2_resid  <- 0.1
shift_Y       <- 2
n_X           <- 1
n_sims        <- as.integer(Sys.getenv("N_SIMS", unset = "500"))
alpha         <- 0.05
n_cores       <- as.integer(Sys.getenv("N_CORES", unset = "10"))
MASTER_SEED   <- 42

tree <- read.tree(TREE_PATH)
phyloKMat <- read.table(KMAT_PATH)
phyloKMat <- phyloKMat / 2
colnames(phyloKMat) <- rownames(phyloKMat)
cat(sprintf("Tree: %d tips (%s)\n", Ntip(tree), TREE_PATH))

# --- helper functions: copied verbatim from src/powerSimulation_XY_revised.R ---

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
  internal_nodes <- internal_nodes[-1]
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
    warning(sprintf("Only found %d independent nodes (requested %d).", length(selected), n))
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

simulate_Y <- function(tree, regime, sigma2_Y, shift_Y, sigma2_resid_Y, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  bm_part    <- fastBM(tree, sig2 = sigma2_Y, internal = FALSE)
  resid_part <- rnorm(Ntip(tree), 0, sqrt(sigma2_resid_Y))
  shift_vec  <- ifelse(regime == 1, shift_Y, 0)
  bm_part + shift_vec + resid_part
}

simulate_X <- function(tree, regime, sigma2_X, shift_X, sigma2_resid_X, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  bm_part    <- fastBM(tree, sig2 = sigma2_X, internal = FALSE)
  resid_part <- rnorm(Ntip(tree), 0, sqrt(sigma2_resid_X))
  shift_vec  <- ifelse(regime == 1, shift_X, 0)
  bm_part + shift_vec + resid_part
}

simulate_Y_and_X <- function(tree, nodes, prop, sigma2_Y, sigma2_X, shift_Y, beta, n_X,
                              sigma2_resid_Y, sigma2_resid_X, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  regime_Y <- make_regime(tree, nodes)
  Y        <- simulate_Y(tree, regime_Y, sigma2_Y, shift_Y, sigma2_resid_Y)

  shift_X  <- beta * shift_Y
  regime_X <- make_regime(tree, sample(nodes, round(length(nodes) * prop)))
  X_mat    <- replicate(n_X, simulate_X(tree, regime_X, sigma2_X, shift_X, sigma2_resid_X))
  colnames(X_mat) <- paste0("X", seq_len(n_X))
  data.frame(assemblyID = as.factor(tree$tip.label), Y = Y, X = X_mat)
}

fit_model <- function(data, Kmat, responseVar, predictor) {
  asreml.options(verbose = FALSE, trace = FALSE, maxit = 60)
  Coeff    <- numeric(length(predictor))
  WaldPval <- numeric(length(predictor))
  for (i in seq_along(predictor)) {
    fullFM     <- as.formula(paste0(responseVar, "~", predictor[i]))
    model_full <- asreml(fixed = fullFM, random = ~ vm(assemblyID, Kmat), ai.sing = FALSE, data = data)
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

run_single_replicate <- function(KMat, tree, nodes, prop, sigma2_Y, sigma2_X, shift_Y, beta, n_X,
                                  sigma2_resid, alpha, sim_seed) {
  sr_Y <- sigma2_resid * sigma2_Y
  sr_X <- sigma2_resid * sigma2_X
  simDat <- simulate_Y_and_X(tree, nodes, prop, sigma2_Y, sigma2_X, shift_Y, beta, n_X, sr_Y, sr_X, seed = sim_seed)
  tmp    <- fit_model(simDat, KMat, "Y", colnames(simDat)[-c(1, 2)])
  hits   <- sum(tmp$p < alpha, na.rm = TRUE)
  total  <- sum(!is.na(tmp$p))
  list(hits = hits, total = total)
}

wilson_ci <- function(hits, n, conf = 0.95) {
  if (n == 0) return(c(lower = NA, upper = NA))
  z     <- qnorm(1 - (1 - conf) / 2)
  p_hat <- hits / n
  denom <- 1 + z^2 / n
  centre <- (p_hat + z^2 / (2 * n)) / denom
  margin <- z * sqrt(p_hat * (1 - p_hat) / n + z^2 / (4 * n^2)) / denom
  c(lower = max(0, centre - margin), upper = min(1, centre + margin))
}

# --- main sweep: one mclapply call over the 15 n_trans values ---

RNGkind("L'Ecuyer-CMRG")
combo_seed <- MASTER_SEED + 999
set.seed(combo_seed)

cat(sprintf("\nSweeping n_trans in {%s} at beta=%.2f, prop=%.0f%%, n_sims=%d\n\n",
            paste(N_TRANS_SEQ, collapse=","), BETA_FIXED, PROP_FIXED*100, n_sims))

results <- mclapply(seq_along(N_TRANS_SEQ), function(ni) {
  n_trans <- N_TRANS_SEQ[ni]
  total_hits  <- 0L
  total_valid <- 0L
  for (s in seq_len(n_sims)) {
    sim_seed  <- combo_seed + ni * 1000 + s
    rep_nodes <- sample_independent_nodes(tree, n_trans, seed = sim_seed)
    res <- tryCatch(
      run_single_replicate(phyloKMat, tree, rep_nodes, PROP_FIXED,
                            sigma2_Y, sigma2_X, shift_Y, BETA_FIXED, n_X,
                            sigma2_resid, alpha, sim_seed),
      error = function(e) {
        warning(sprintf("Replicate failed (n_trans=%d, sim=%d): %s", n_trans, s, conditionMessage(e)))
        list(hits = 0L, total = 0L)
      }
    )
    total_hits  <- total_hits  + res$hits
    total_valid <- total_valid + res$total
  }
  power <- total_hits / total_valid
  ci    <- wilson_ci(total_hits, total_valid)
  cat(sprintf("n_trans=%2d: power=%.3f [%.3f, %.3f]  (%d/%d)\n",
              n_trans, power, ci[1], ci[2], total_hits, total_valid))
  c(n_trans = n_trans, power = power, ci_lower = ci[1], ci_upper = ci[2],
    hits = total_hits, total = total_valid)
}, mc.cores = n_cores)

df <- as.data.frame(do.call(rbind, results))
df$beta <- BETA_FIXED
df$prop <- PROP_FIXED

write.csv(df, file.path(OUTPUT_DIR, "power_vs_transitionCount.csv"), row.names = FALSE)

cat("\n--- Power vs. n_trans (beta=0.3, 100% overlap) ---\n")
print(df, row.names = FALSE)

crossing <- df[df$power >= POWER_TARGET, ]
if (nrow(crossing) > 0) {
  n_needed <- min(crossing$n_trans)
  cat(sprintf("\nSmallest n_trans reaching power >= %.2f: %d (power = %.3f)\n",
              POWER_TARGET, n_needed, df$power[df$n_trans == n_needed]))
} else {
  cat(sprintf("\nPower target %.2f not reached within the tested range (max n_trans=%d, power=%.3f)\n",
              POWER_TARGET, max(df$n_trans), max(df$power)))
}

cat(sprintf("\nWrote %s\n", file.path(OUTPUT_DIR, "power_vs_transitionCount.csv")))
