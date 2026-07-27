# =============================================================================
# Mediation Analysis: Disease (Delirium) → Cholinergic Dysfunction → Segregation
# =============================================================================
# Model: X (delirium) → M (cholinergic variable) → Y (mean participation coefficient)
#
# Required packages: mediation, dplyr, ggplot2
# Install if needed: install.packages(c("mediation", "dplyr", "ggplot2"))
# =============================================================================

library(mediation)
library(dplyr)
library(ggplot2)

# ---- 1. Load Data -----------------------------------------------------------
setwd('/Users/ntaylor/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)/Postdoc_Rob/Analysis/Graph_Theory/schaef_400/csf_regress_redo')
df <- read.csv("mediation_variables_postpeak_pc.csv", stringsAsFactors = FALSE)
df$delirium <- as.factor(df$delirium)

# ── Define combinations ──
seeds <- c("LC", "nbM")
mediators <- c("peak_freq", "std_seed", "mean_peak_amp")
mediator_labels <- c("Peak Frequency", "Seed SD", "Mean Peak Amplitude")

# ── Storage ──
results_list <- list()

for (seed in seeds) {
  df_seed <- df %>% filter(seed_name == seed)
  
  cat(paste0("\n===========================================================\n"))
  cat(paste0("  SEED: ", seed, "  (n = ", nrow(df_seed), ")\n"))
  cat(paste0("===========================================================\n"))
  
  for (m in seq_along(mediators)) {
    med_var <- mediators[m]
    med_label <- mediator_labels[m]
    
    cat(paste0("\n--- Mediator: ", med_label, " (", med_var, ") ---\n"))
    
    # ── Path c: Total effect (Disease → Post-Peak PC) ──
    model_total <- lm(
      as.formula("mean_pc_post_peak ~ delirium"),
      data = df_seed
    )
    cat("\n  Path c (Total Effect: Delirium → Post-Peak PC):\n")
    print(summary(model_total)$coefficients)
    
    # ── Path a: Disease → Mediator ──
    model_m <- lm(
      as.formula(paste0(med_var, " ~ delirium")),
      data = df_seed
    )
    cat(paste0("\n  Path a (Delirium → ", med_label, "):\n"))
    print(summary(model_m)$coefficients)
    
    # ── Paths b & c': Disease + Mediator → Post-Peak PC ──
    model_y <- lm(
      as.formula(paste0("mean_pc_post_peak ~ delirium + ", med_var)),
      data = df_seed
    )
    cat(paste0("\n  Paths b & c' (Delirium + ", med_label, " → Post-Peak PC):\n"))
    print(summary(model_y)$coefficients)
    
    # ── Mediation ──
    cat("Running nonparametric bootstrap\n")
    med_result <- mediate(
      model_m, model_y,
      treat = "delirium",
      mediator = med_var,
      treat.value = "1",
      control.value = "0",
      boot = TRUE,
      sims = 5000,
      boot.ci.type = "perc"
    )
    
    cat(paste0("\n  Mediation Results (", med_label, 
               ", 5000 bootstrap iterations):\n"))
    print(summary(med_result))
    
    # ── Store results ──
    s <- summary(med_result)
    results_list[[paste0(seed, "_", med_var)]] <- data.frame(
      seed = seed,
      mediator = med_var,
      mediator_label = med_label,
      outcome = "mean_pc_post_peak",
      acme_estimate = s$d0,
      acme_ci_lower = s$d0.ci[1],
      acme_ci_upper = s$d0.ci[2],
      acme_p = s$d0.p,
      ade_estimate = s$z0,
      ade_ci_lower = s$z0.ci[1],
      ade_ci_upper = s$z0.ci[2],
      ade_p = s$z0.p,
      total_estimate = s$tau.coef,
      total_ci_lower = s$tau.ci[1],
      total_ci_upper = s$tau.ci[2],
      total_p = s$tau.p,
      prop_mediated = s$n0,
      prop_med_ci_lower = s$n0.ci[1],
      prop_med_ci_upper = s$n0.ci[2],
      prop_med_p = s$n0.p,
      stringsAsFactors = FALSE
    )
  }
}

# ── Summary table ──
summary_df <- bind_rows(results_list)
write.csv(summary_df, "mediation_results_postpeak.csv", row.names = FALSE)
cat("\nSaved: mediation_results_postpeak.csv\n")

# ── Forest plot ──
plot_df <- summary_df %>%
  mutate(
    model_label = reorder(paste0(seed, "\n", mediator_label), acme_estimate),
    sig = ifelse(acme_p < 0.05, "Significant", "Non-significant")
  )

p <- ggplot(plot_df, aes(x = acme_estimate, y = model_label)) +
  geom_point(aes(colour = sig), size = 3) +
  geom_errorbarh(aes(xmin = acme_ci_lower, xmax = acme_ci_upper, colour = sig),
                 height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "red") +
  scale_colour_manual(values = c("Significant" = "black",
                                  "Non-significant" = "grey60"),
                      name = "") +
  labs(
    title = "Indirect Effects (ACME) — Post-Peak Participation Coefficient",
    subtitle = "Disease → Cholinergic Mediator → Post-Peak Network Segregation",
    x = "Indirect Effect (ACME) with 95% Bootstrap CI",
    y = ""
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 11),
    legend.position = "bottom"
  )

ggsave("forest_plot_mediation_postpeak.png", p,
       width = 8, height = 5, dpi = 300)
cat("Saved: forest_plot_mediation_postpeak.png\n")

# ── Sensitivity analysis for significant results ──
for (key in names(results_list)) {
  r <- results_list[[key]]
  if (r$acme_p < 0.05) {
    cat(paste0("\n=== Sensitivity Analysis: ", key, " ===\n"))
    
    df_seed <- df %>% filter(seed_name == r$seed)
    
    model_m <- lm(
      as.formula(paste0(r$mediator, " ~ delirium")),
      data = df_seed
    )
    model_y <- lm(
      as.formula(paste0("mean_pc_post_peak ~ delirium + ", r$mediator)),
      data = df_seed
    )
    
    med_result <- mediate(
      model_m, model_y,
      treat = "delirium",
      mediator = r$mediator,
      treat.value = "1",
      control.value = "0",
      boot = TRUE,
      sims = 5000,
      boot.ci.type = "perc"
    )
    
    sens <- medsens(med_result, rho.by = 0.01)
    print(summary(sens))
    
    png(paste0("sensitivity_", key, ".png"), width = 800, height = 600)
    plot(sens, main = paste0("Sensitivity: ", r$seed, " - ", r$mediator_label))
    dev.off()
    cat(paste0("Saved: sensitivity_", key, ".png\n"))
  }
}