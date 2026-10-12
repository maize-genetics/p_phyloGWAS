# Power curve for the empirical-transition scenario, one line per convergence level
# (redraw of Fig2_empiricalTransition_byConvergence, without the title).
#
# Input: the 687-tip power simulation (src/powerSimulation_XY_revised.R run per scenario
# into output/revisionAnalyses/powerSimulation_687/<scenario>/power_df_merged.csv; see
# power_data() below) and output/powerSimulation/realTransition.json (687-tip tree).
# Defines power_curve_plot(), which 12M_powerSimIllustration.R sources to build the
# combined figure (Fig 2B-C), and power_curve_by_scenario() for the supplementary
# figure. Run directly, it saves the standalone curve and the supplementary figure.
#
# Usage: Rscript notebook/12_revisionAnalyses/12M_powerCurve.R

suppressMessages({library(ggplot2)})

PHYLOGWAS_ROOT <- Sys.getenv("PHYLOGWAS_ROOT", unset = "/workdir/sh2246/p_phyloGWAS")

# Power results on the 687-tip tree: one power_df_merged.csv per scenario job.
POWER_DIR <- file.path(PHYLOGWAS_ROOT, "output/revisionAnalyses/powerSimulation_687")
power_data <- function(scenarios = NULL) {
  f <- list.files(POWER_DIR, pattern = "^power_df_merged\\.csv$", recursive = TRUE, full.names = TRUE)
  d <- do.call(rbind, lapply(f, read.csv))
  if (!is.null(scenarios)) {
    miss <- setdiff(scenarios, d$scenario)
    if (length(miss)) stop("power results missing for: ", paste(miss, collapse = ", "))
    d <- d[d$scenario %in% scenarios, ]
  }
  d
}

# Sequential (viridis) colours: convergence is ordered, so lightness encodes it.
CONV_PAL <- c("100%" = "#440154", "80%" = "#3B528B", "50%" = "#21918C", "30%" = "#5EC962")

# base_size in pt. legend_inside puts the legend in the empty lower-right corner,
# which saves width when the panel sits next to another one.
power_curve_plot <- function(base_size = 13, legend_inside = FALSE, title = NULL) {
  e <- power_data("empirical transition")
  e$convergence <- factor(e$convergence, levels = names(CONV_PAL))
  k <- base_size / 13                       # scale line and point sizes with the text
  # Legend: convergence level and the number of co-shift branches it implies, i.e.
  # round(n_transitions * level), the same rounding the simulation uses.
  n_trans <- length(unlist(rjson::fromJSON(
    file = file.path(PHYLOGWAS_ROOT, "output/powerSimulation/realTransition.json"))[[1]]))
  conv_lab <- sprintf("%s (%d)", names(CONV_PAL),
                      round(n_trans * as.numeric(sub("%", "", names(CONV_PAL))) / 100))
  ggplot(e, aes(beta, power, colour = convergence, fill = convergence)) +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = .3, colour = NA) +
    geom_hline(yintercept = .8, linetype = "dotted", colour = "black", linewidth = .8 * k) +
    geom_line(linewidth = 1.2 * k) +
    geom_point(size = 2.5 * k) +
    scale_colour_manual(values = CONV_PAL, labels = conv_lab) +
    scale_fill_manual(values = CONV_PAL, labels = conv_lab) +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    labs(x = expression(theta ~ "magnitude"), y = "Estimated power",
         colour = "Convergence\n(co-shift branches)", fill = "Convergence\n(co-shift branches)",
         title = title) +
    theme_classic(base_size = base_size) +
    theme(legend.position = if (legend_inside) c(.98, .03) else "right",
          legend.justification = if (legend_inside) c(1, 0) else "center",
          legend.background = element_blank(),
          plot.title = element_text(size = base_size, hjust = 0, margin = margin(b = 5.5, l = 8)),
          plot.title.position = "plot",
          legend.key.size = grid::unit(base_size * .9, "pt"))
}

# Supplementary version: all four transition scenarios as facets, convergence levels
# as lines (same style as power_curve_plot). Each facet is labelled with the number
# of transition branches its scenario draws; co-shift branches = that x convergence.
power_curve_by_scenario <- function(base_size = 8) {
  d <- power_data(c("empirical transition", "random transition - 1/2x",
                    "random transition - 1x", "random transition - 2x"))
  d$convergence <- factor(d$convergence, levels = names(CONV_PAL))
  n_trans <- length(unlist(rjson::fromJSON(
    file = file.path(PHYLOGWAS_ROOT, "output/powerSimulation/realTransition.json"))[[1]]))
  sc <- c("empirical transition" = 1, "random transition - 1/2x" = .5,
          "random transition - 1x" = 1, "random transition - 2x" = 2)
  sc_lab <- c(sprintf("Empirical transitions (%d branches)", n_trans),
              sprintf("Random, \u00bd\u00d7 (%d branches)", round(n_trans * .5)),
              sprintf("Random, 1\u00d7 (%d branches)", n_trans),
              sprintf("Random, 2\u00d7 (%d branches)", n_trans * 2))
  # round(n_trans * .5) matches the simulation's own rounding for the 1/2x scenario.
  d$scenario <- factor(d$scenario, levels = names(sc), labels = sc_lab)
  k <- base_size / 13
  ggplot(d, aes(beta, power, colour = convergence, fill = convergence)) +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = .3, colour = NA) +
    geom_hline(yintercept = .8, linetype = "dotted", colour = "black", linewidth = .8 * k) +
    geom_line(linewidth = 1.2 * k) +
    geom_point(size = 2.5 * k) +
    facet_wrap(~ scenario, nrow = 2) +
    scale_colour_manual(values = CONV_PAL) + scale_fill_manual(values = CONV_PAL) +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    labs(x = expression(theta ~ "magnitude"), y = "Estimated power",
         colour = "Convergence", fill = "Convergence") +
    theme_classic(base_size = base_size) +
    theme(strip.background = element_blank(), strip.text = element_text(size = base_size, hjust = 0),
          legend.background = element_blank(), panel.spacing = grid::unit(10, "pt"),
          legend.key.size = grid::unit(base_size * .9, "pt"))
}

if (sys.nframe() == 0) {
  out <- file.path(PHYLOGWAS_ROOT, "output/revisionAnalyses/powerSimulation/powerCurve_empirical_byConvergence.png")
  dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
  ggsave(out, power_curve_plot(), width = 8, height = 5.5, dpi = 300)
  cat("wrote", out, "\n")
  # Supplementary figure (replaces the convergence-faceted Fig2_revised layout).
  stem <- file.path(dirname(out), "FigS_powerCurve_byScenario_17.4cm")
  # Needs all four scenario runs; skipped (with a message) until they are done.
  fs <- tryCatch(power_curve_by_scenario(), error = function(e) { message("skip supplementary figure: ", conditionMessage(e)); NULL })
  if (!is.null(fs)) {
    ggsave(paste0(stem, ".png"), fs, width = 17.4, height = 11, units = "cm", dpi = 600)
    ggsave(paste0(stem, ".pdf"), fs, width = 17.4, height = 11, units = "cm", device = cairo_pdf)
    cat("wrote", stem, "{png,pdf}\n")
  }
}
