# Illustration of the power-simulation model (src/powerSimulation_XY_revised.R) on a
# toy tree (50 tips, 10 independent transition branches), combined with the empirical
# power curve (12M_powerCurve.R) into one figure.
#
# Same generative model as the simulation:
#   Y = BM(sigma2_Y) + shift_Y on transition clades + residual
#   X = BM(sigma2_X) + theta * shift_Y on a fraction (convergence) of those clades + residual
#
# Columns: Y; X with no co-shift (theta = 0); X with theta = 0.3 and 0.5 at 100%
# convergence; X with theta = 0.5 at 30% convergence. The tree, the transition branches and every
# BM / residual draw are held fixed across columns, so the only thing that changes
# between columns is the X shift. Shifted tips also show their neutral expectation
# (the value without the shift) in light grey.
#
# Outputs (output/revisionAnalyses/powerSimulation/), all drawn natively at print size:
#   powerSim_illustration_50tips_10transitions.png       illustration alone
#   Fig_powerSim_horizontal_17.4cm.{png,pdf}              A | B, 17.4 x 8.7 cm
#   Fig_powerSim_vertical_17.4cm.{png,pdf}                A / B, 17.4 x 13 cm
#
# Usage: Rscript notebook/12_revisionAnalyses/12M_powerSimIllustration.R

suppressMessages({
  library(ape); library(phytools); library(ggtree); library(ggplot2)
  library(patchwork); library(data.table)
})

PHYLOGWAS_ROOT <- Sys.getenv("PHYLOGWAS_ROOT", unset = "/workdir/sh2246/p_phyloGWAS")
OUT_DIR <- file.path(PHYLOGWAS_ROOT, "output/revisionAnalyses/powerSimulation")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
source(file.path(PHYLOGWAS_ROOT, "notebook/12_revisionAnalyses/12M_powerCurve.R"))

# ---- parameters (simulation values) ----
N_TIP <- 50; N_TRANS <- 10; SEED <- as.integer(Sys.getenv("SEED", unset = "7"))
sigma2_Y <- 0.5; sigma2_X <- 0.5; sigma2_resid <- 0.1; shift_Y <- 2
X_SETTINGS <- list(c(0, 1), c(0.3, 1), c(0.5, 1), c(0.5, 0.3))   # X columns as (theta, convergence)

# ---- helpers, as in src/powerSimulation_XY_revised.R ----
get_all_ancestors <- function(tree, node) {
  anc <- c(); cur <- node; root <- Ntip(tree) + 1
  while (cur != root) { par <- tree$edge[tree$edge[, 2] == cur, 1]; if (!length(par)) break
                        anc <- c(anc, par); cur <- par }
  anc
}
sample_independent_nodes <- function(tree, n) {
  internal <- (Ntip(tree) + 2):(Ntip(tree) + Nnode(tree))   # exclude root
  sel <- c()
  for (nd in sample(internal)) {
    if (length(sel) == n) break
    if (!any(sel %in% phytools::getDescendants(tree, nd)) && !any(sel %in% get_all_ancestors(tree, nd)))
      sel <- c(sel, nd)
  }
  sel
}
tips_below <- function(tree, nd) { d <- phytools::getDescendants(tree, nd); tree$tip.label[d[d <= Ntip(tree)]] }
regime_of <- function(tree, nodes) {
  r <- setNames(rep(0L, Ntip(tree)), tree$tip.label)
  for (nd in nodes) r[tips_below(tree, nd)] <- 1L
  r
}

# ---- one fixed draw of everything ----
set.seed(SEED)
tree <- rcoal(N_TIP, tip.label = sprintf("t%02d", 1:N_TIP))
tree$edge.length <- tree$edge.length / max(node.depth.edgelength(tree))   # height 1
nodes <- sample_independent_nodes(tree, N_TRANS)
conv_order <- sample(nodes)               # X joins clades in this order as convergence rises
regY <- regime_of(tree, nodes)

bmY <- fastBM(tree, sig2 = sigma2_Y); eY <- rnorm(N_TIP, 0, sqrt(sigma2_resid))
bmX <- fastBM(tree, sig2 = sigma2_X); eX <- rnorm(N_TIP, 0, sqrt(sigma2_resid))
Y <- bmY + shift_Y * regY + eY
make_X <- function(theta, prop) {
  regX <- regime_of(tree, conv_order[seq_len(round(N_TRANS * prop))])
  list(value = bmX + theta * shift_Y * regX + eX, shifted = regX * (theta > 0))
}
Xs <- lapply(X_SETTINGS, function(st) make_X(st[1], st[2]))
titles <- sapply(X_SETTINGS, function(st)
  if (st[1] == 0) "X, θ = 0\nno co-shift" else sprintf("X, θ = %.2g\n%d%% convergent", st[1], round(100 * st[2])))

# ---- drawing; every size scales with base font size fs (pt) ----
# Unshifted tips are mid grey. A shifted tip is drawn at its shifted value (colour)
# with its neutral expectation (the same tip without the shift) in light grey, joined
# by a coloured segment, so the shift reads as the gap between the two.
trans_col <- "#c0392b"; conv_col <- "#2471a3"; base_col <- "grey55"; ghost_col <- "grey82"
states <- c("background", "Y transition", "X co-shift")
pal <- c(background = base_col, `Y transition` = trans_col, `X co-shift` = conv_col)
ylims <- c(.5, N_TIP + .5)
xlims <- range(c(Y, bmX + eX, bmX + eX + shift_Y * max(sapply(X_SETTINGS, `[`, 1)))) + c(-.2, .2)
tip_y <- {p0 <- ggtree(tree, branch.length = "none"); d0 <- p0$data[p0$data$isTip, ]; setNames(d0$y, d0$label)}

illustration <- function(fs, title = NULL, title_size = fs) {
  k <- fs / 9
  ttl <- element_text(size = fs * .9, hjust = .5, lineheight = .95)
  # Cladogram, so short transition branches near the tips stay visible.
  tr <- {
    p <- ggtree(tree, colour = "grey35", linewidth = .35 * k, branch.length = "none")
    d <- p$data
    d$hl <- d$node %in% unlist(lapply(nodes, function(nd) c(nd, phytools::getDescendants(tree, nd))))
    p + geom_tree(data = d[d$hl, ], colour = trans_col, linewidth = .7 * k) +
      geom_point2(aes(subset = node %in% nodes), colour = trans_col, size = 2 * k, shape = 18) +
      scale_y_continuous(limits = ylims, expand = c(0, 0)) +
      labs(title = " ") + theme(plot.title = ttl)   # blank title keeps the panel aligned
  }
  strip <- function(value, neutral, state, title) {
    dt <- data.table(value = value, neutral = neutral[names(value)], y = tip_y[names(value)],
                     state = factor(state, states))
    sh <- dt[state != "background"]
    ggplot(dt, aes(value, y, colour = state)) +
      geom_vline(xintercept = 0, colour = "grey85", linewidth = .3 * k) +
      geom_segment(aes(x = 0, xend = neutral, yend = y), colour = "grey88", linewidth = .3 * k) +
      geom_point(data = sh, aes(x = neutral), colour = ghost_col, size = 1.1 * k) +
      geom_segment(data = sh, aes(x = neutral, xend = value, yend = y), linewidth = .45 * k) +
      geom_point(size = 1.1 * k) +
      scale_colour_manual(values = pal, drop = FALSE) +
      scale_y_continuous(limits = ylims, expand = c(0, 0), breaks = NULL) +
      scale_x_continuous(limits = xlims, breaks = c(-2, 0, 2)) +
      labs(x = NULL, y = NULL, title = title) +
      theme_bw(fs) + theme(plot.title = ttl, panel.grid = element_blank(), legend.position = "none",
                           axis.text = element_text(size = fs * .8),
                           # transparent, so a title wider than its column is not painted over
                           plot.background = element_blank())
  }
  # Colour key, as a narrow column at the right.
  key_dt <- data.frame(s = factor(c(states, "background"), states), x = 0, y = c(.68, .56, .44, .32),
                       col = c(base_col, trans_col, conv_col, ghost_col),
                       label = c("no shift", "Y transition", "X co-shift", "neutral\nexpectation"))
  key <- ggplot(key_dt, aes(x, y)) +
    geom_point(aes(colour = I(col)), size = 1.6 * k) +
    geom_text(aes(x = .12, label = label), hjust = 0, size = fs * .85 / .pt, lineheight = .9, colour = "black") + coord_cartesian(xlim = c(-.05, 1), ylim = c(0, 1), clip = "off") +
    theme_void() + theme(legend.position = "none")

  cells <- c(list(T = tr, Y = strip(Y, bmY + eY, ifelse(regY == 1, "Y transition", "background"), "Y\n(environment)")),
             setNames(Map(function(X, t) strip(X$value, bmX + eX, ifelse(X$shifted == 1, "X co-shift", "background"), t),
                          Xs, titles), letters[seq_along(Xs)]),
             list(L = key))
  wrap_plots(cells, design = paste0("TY", paste(letters[seq_along(Xs)], collapse = ""), "L"),
             widths = c(1.15, 1, rep(1, length(Xs)), .9)) +
    plot_annotation(title = title,
                    theme = theme(plot.title = element_text(size = title_size, hjust = 0,
                                                            margin = margin(b = 5.5, l = 8)),
                                  plot.margin = margin(5.5, 5.5, 5.5, 0)))   # title flush with panel B's
}

# ---- illustration alone ----
out <- file.path(OUT_DIR, "powerSim_illustration_50tips_10transitions.png")
ggsave(out, illustration(10), width = 11, height = 6, dpi = 300)
cat("wrote", out, "\n")

# ---- combined figures at journal width (17.4 cm) ----
FS <- 7   # base font size, pt
panel <- function(p) wrap_elements(full = p)
save_both <- function(fig, stem, w_cm, h_cm) {
  for (ext in c("png", "pdf")) {
    f <- file.path(OUT_DIR, sprintf("%s.%s", stem, ext))
    if (ext == "pdf") ggsave(f, fig, width = w_cm, height = h_cm, units = "cm", device = cairo_pdf)
    else ggsave(f, fig, width = w_cm, height = h_cm, units = "cm", dpi = 600)
  }
  cat("wrote", file.path(OUT_DIR, stem), sprintf("{png,pdf}  %.1f x %.1f cm\n", w_cm, h_cm))
}

# Panel B text is one size step (1 pt) larger than panel A.
B_TITLE <- "Power analysis based on empirical phylogeny and transitions"
A <- panel(illustration(FS, title = "Simplified illustration of the simulation", title_size = FS + 1))
# Horizontal: legend inside the curve keeps panel B narrow.
B_h <- panel(power_curve_plot(FS + 1, legend_inside = TRUE, title = B_TITLE))
save_both(A + B_h + plot_layout(widths = c(2.7, 1)), "Fig_powerSim_horizontal_17.4cm", 17.4, 8.7)
# Vertical: panel B gets the full width, legend at right.
B_v <- panel(power_curve_plot(FS + 1, title = B_TITLE))
save_both(A / B_v + plot_layout(heights = c(1.1, 1)), "Fig_powerSim_vertical_17.4cm", 17.4, 13)
