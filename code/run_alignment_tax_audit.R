#!/usr/bin/env Rscript

# Minimal end-to-end audit for the archived manuscript.
# The script:
# 1. filters the April 8, 2024 item-level leaderboard snapshot,
# 2. computes the PC1-based ability metric,
# 3. merges precomputed jackknife influence values,
# 4. fits the reported regression models, and
# 5. writes the figures and tables consumed by the manuscript.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
})

`%||%` <- function(a, b) {
  if (!is.null(a)) a else b
}

parse_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (!startsWith(arg, "--")) next
    kv <- sub("^--", "", arg)
    if (grepl("=", kv, fixed = TRUE)) {
      parts <- strsplit(kv, "=", fixed = TRUE)[[1]]
      out[[parts[1]]] <- parts[2]
    } else {
      out[[kv]] <- TRUE
    }
  }
  out
}

cosine_similarity <- function(a, b) {
  denom <- sqrt(sum(a * a)) * sqrt(sum(b * b))
  if (denom == 0) return(NA_real_)
  sum(a * b) / denom
}

classify_model_type <- function(model_id) {
  id <- tolower(model_id)
  if (str_detect(id, "merge|mixture|mix|blend|passthrough|fusion|fuse|slerp|grafted|ties")) {
    return("merge_or_moe")
  }
  if (str_detect(id, "instruct|chat|alpaca|orca|tulu|rlhf|dpo|guanaco|assistant|openbuddy|it\\b")) {
    return("instruction_or_chat")
  }
  "base_or_other"
}

coef_table_with_ci <- function(fit) {
  sm <- summary(fit)$coefficients
  est <- sm[, 1]
  se <- sm[, 2]
  t_crit <- qt(0.975, df = fit$df.residual)
  tibble(
    term = rownames(sm),
    estimate = est,
    std_error = se,
    t_value = sm[, 3],
    p_value = sm[, 4],
    ci_low_95 = est - t_crit * se,
    ci_high_95 = est + t_crit * se
  )
}

safe_fit_r2 <- function(formula, data) {
  fit <- tryCatch(lm(formula, data = data), error = function(e) NULL)
  if (is.null(fit)) return(NA_real_)
  r2 <- tryCatch(summary(fit)$r.squared, error = function(e) NA_real_)
  if (!is.finite(r2)) return(NA_real_)
  r2
}

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else ""
script_dir <- if (nzchar(script_path)) dirname(normalizePath(script_path)) else getwd()
project_root <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)
args <- parse_args(commandArgs(trailingOnly = TRUE))

input_csv <- args$input %||% file.path(project_root, "data", "hugging_item_level_data_final_april8_2024.csv")
filter_csv <- args$filter_csv %||% file.path(project_root, "data", "test_level_data_july23_2024_Levenshtein_20_N_641_factors.csv")
jackknife_csv <- args$jackknife_csv %||% file.path(project_root, "data", "jackknife_results_checkpoint.csv")
outdir <- args$outdir %||% file.path(project_root, "outputs")
recompute_jackknife <- isTRUE(args$recompute_jackknife) || identical(args$recompute_jackknife, "true")
bootstrap_b <- as.integer(args$bootstrap_b %||% 10000)
seed <- as.integer(args$seed %||% 20260223)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
set.seed(seed)

cat("== Alignment Tax Audit ==\n")
cat("input_csv:", input_csv, "\n")
cat("filter_csv:", filter_csv, "\n")
cat("jackknife_csv:", jackknife_csv, "\n")
cat("outdir:", outdir, "\n")
cat("recompute_jackknife:", recompute_jackknife, "\n")
cat("bootstrap_b:", bootstrap_b, "\n")

raw <- read_csv(input_csv, show_col_types = FALSE)
flt <- read_csv(filter_csv, show_col_types = FALSE)

id_col <- if ("model_id" %in% names(flt)) "model_id" else "Model"
ids <- flt[[id_col]]

dat <- raw %>%
  filter(model_id %in% ids) %>%
  distinct(model_id, .keep_all = TRUE)

if (nrow(dat) == 0) {
  stop("No rows after filtering. Check input/filter paths and ID columns.")
}

model_ids <- dat$model_id
X <- as.matrix(dat %>% select(-model_id))
storage.mode(X) <- "numeric"

N <- nrow(X)
P <- ncol(X)
cat("filtered_models:", N, "\n")
cat("items:", P, "\n")

# Item-centered PCA. `prcomp` handles full SVD with deterministic orientation up to sign.
pca <- prcomp(X, center = TRUE, scale. = FALSE)
eig <- pca$sdev^2
var_pct <- eig / sum(eig)

meta <- flt %>%
  transmute(
    model_id = .data[[id_col]],
    params_num = if ("params_num" %in% names(flt)) params_num else NA_real_,
    gfactor_meta = if ("gfactor" %in% names(flt)) gfactor else NA_real_
  )

tmp <- tibble(model_id = model_ids, ability_raw = pca$x[, 1]) %>%
  left_join(meta, by = "model_id")

sign_flip <- ifelse(cor(tmp$ability_raw, tmp$gfactor_meta, use = "pairwise.complete.obs") < 0, -1, 1)
if (is.na(sign_flip)) sign_flip <- -1

ability_pc1 <- sign_flip * pca$x[, 1]
ability <- 100 + 15 * (ability_pc1 - mean(ability_pc1)) / sd(ability_pc1)
loading1 <- sign_flip * pca$rotation[, 1]

mu <- colMeans(X)
X_centered <- sweep(X, 2, mu, "-")

coh_cos_raw <- apply(X, 1, function(row) cosine_similarity(row, loading1))
coh_cos_centered <- apply(X_centered, 1, function(row) cosine_similarity(row, loading1))

Xhat1 <- tcrossprod(ability_pc1, loading1) + matrix(mu, nrow = N, ncol = P, byrow = TRUE)
sq_resid <- (X - Xhat1)^2
rmse_1f <- sqrt(rowMeans(sq_resid))
ss_total <- rowSums((X - matrix(mu, nrow = N, ncol = P, byrow = TRUE))^2)
r2_1f <- 1 - rowSums(sq_resid) / ss_total

baseline_pc1_var <- var_pct[1]

if (recompute_jackknife) {
  # Recomputing all leave-one-out PCA fits is slow, so the archived workflow
  # defaults to the checkpointed jackknife results used for the paper.
  cat("Computing full jackknife from scratch (this can take ~1 hour on this hardware)...\n")
  delta <- numeric(N)
  for (i in seq_len(N)) {
    pca_i <- prcomp(X[-i, , drop = FALSE], center = TRUE, scale. = FALSE, rank. = 1)
    var_i <- (pca_i$sdev[1]^2) / sum(pca_i$sdev^2)
    delta[i] <- var_i - baseline_pc1_var
    if (i %% 25 == 0 || i == N) {
      cat(sprintf("  jackknife %d/%d\n", i, N))
    }
  }
  jackknife <- tibble(Model = model_ids, Delta_g = delta)
  write_csv(jackknife, file.path(outdir, "jackknife_recomputed.csv"))
} else {
  if (!file.exists(jackknife_csv)) {
    stop("jackknife_csv not found; rerun with --recompute_jackknife=true or provide a valid file.")
  }
  jackknife <- read_csv(jackknife_csv, show_col_types = FALSE) %>%
    select(Model, Delta_g)
}

metrics <- tibble(
  Model = model_ids,
  Ability_PC1 = ability_pc1,
  Ability = ability,
  Coherence_Cosine_Raw = coh_cos_raw,
  Coherence_Cosine_Centered = coh_cos_centered,
  Coherence_R2_1F = r2_1f,
  RMSE_1F = rmse_1f
) %>%
  left_join(jackknife, by = "Model") %>%
  left_join(meta, by = c("Model" = "model_id")) %>%
  mutate(
    model_type = vapply(Model, classify_model_type, FUN.VALUE = character(1)),
    size_band = case_when(
      is.na(params_num) ~ "unknown",
      params_num < 2 ~ "<2B",
      params_num < 20 ~ "2B-20B",
      TRUE ~ "20B+"
    ),
    ability_decile = ntile(Ability, 10)
  )

# Main manuscript models.
linear_fit <- lm(Delta_g ~ Ability, data = metrics)
quad_fit <- lm(Delta_g ~ Ability + I(Ability^2), data = metrics)
full_fit <- lm(Delta_g ~ Ability + I(Ability^2) + log1p(params_num) + model_type, data = metrics)

anova_lq <- anova(linear_fit, quad_fit)
linear_coef_tbl <- coef_table_with_ci(linear_fit)
quad_coef_tbl <- coef_table_with_ci(quad_fit)
full_coef_tbl <- coef_table_with_ci(full_fit)

quad_beta <- coef(quad_fit)
turning_point_est <- -quad_beta[["Ability"]] / (2 * quad_beta[["I(Ability^2)"]])
peak_delta_est <- quad_beta[["(Intercept)"]] +
  quad_beta[["Ability"]] * turning_point_est +
  quad_beta[["I(Ability^2)"]] * turning_point_est^2

boot_stats <- vector("list", bootstrap_b)
for (b in seq_len(bootstrap_b)) {
  idx <- sample.int(N, N, replace = TRUE)
  d <- metrics[idx, ]
  linear_fit_b <- lm(Delta_g ~ Ability, data = d)
  quad_fit_b <- lm(Delta_g ~ Ability + I(Ability^2), data = d)
  qb <- coef(quad_fit_b)
  tp <- if (!is.na(qb[["I(Ability^2)"]]) && qb[["I(Ability^2)"]] != 0) {
    -qb[["Ability"]] / (2 * qb[["I(Ability^2)"]])
  } else {
    NA_real_
  }
  peak <- if (!is.na(tp)) {
    qb[["(Intercept)"]] + qb[["Ability"]] * tp + qb[["I(Ability^2)"]] * tp^2
  } else {
    NA_real_
  }
  boot_stats[[b]] <- tibble(
    b = b,
    cor_ability_delta = cor(d$Ability, d$Delta_g, use = "pairwise.complete.obs"),
    cor_ability_r2 = cor(d$Ability, d$Coherence_R2_1F, use = "pairwise.complete.obs"),
    cor_delta_r2 = cor(d$Delta_g, d$Coherence_R2_1F, use = "pairwise.complete.obs"),
    linear_r2 = summary(linear_fit_b)$r.squared,
    quad_r2 = summary(quad_fit_b)$r.squared,
    full_r2 = safe_fit_r2(Delta_g ~ Ability + I(Ability^2) + log1p(params_num) + model_type, d),
    turning_point = tp,
    peak_delta = peak
  )
}
boot_tbl <- bind_rows(boot_stats)

boot_ci <- boot_tbl %>%
  summarise(
    cor_ability_delta_mean = mean(cor_ability_delta, na.rm = TRUE),
    cor_ability_delta_lo = quantile(cor_ability_delta, 0.025, na.rm = TRUE),
    cor_ability_delta_hi = quantile(cor_ability_delta, 0.975, na.rm = TRUE),
    cor_ability_r2_mean = mean(cor_ability_r2, na.rm = TRUE),
    cor_ability_r2_lo = quantile(cor_ability_r2, 0.025, na.rm = TRUE),
    cor_ability_r2_hi = quantile(cor_ability_r2, 0.975, na.rm = TRUE),
    cor_delta_r2_mean = mean(cor_delta_r2, na.rm = TRUE),
    cor_delta_r2_lo = quantile(cor_delta_r2, 0.025, na.rm = TRUE),
    cor_delta_r2_hi = quantile(cor_delta_r2, 0.975, na.rm = TRUE),
    linear_r2_mean = mean(linear_r2, na.rm = TRUE),
    linear_r2_lo = quantile(linear_r2, 0.025, na.rm = TRUE),
    linear_r2_hi = quantile(linear_r2, 0.975, na.rm = TRUE),
    quad_r2_mean = mean(quad_r2, na.rm = TRUE),
    quad_r2_lo = quantile(quad_r2, 0.025, na.rm = TRUE),
    quad_r2_hi = quantile(quad_r2, 0.975, na.rm = TRUE),
    full_r2_mean = mean(full_r2, na.rm = TRUE),
    full_r2_lo = quantile(full_r2, 0.025, na.rm = TRUE),
    full_r2_hi = quantile(full_r2, 0.975, na.rm = TRUE),
    turning_point_mean = mean(turning_point, na.rm = TRUE),
    turning_point_lo = quantile(turning_point, 0.025, na.rm = TRUE),
    turning_point_hi = quantile(turning_point, 0.975, na.rm = TRUE),
    peak_delta_mean = mean(peak_delta, na.rm = TRUE),
    peak_delta_lo = quantile(peak_delta, 0.025, na.rm = TRUE),
    peak_delta_hi = quantile(peak_delta, 0.975, na.rm = TRUE)
  )

decile_tbl <- metrics %>%
  group_by(ability_decile) %>%
  summarise(
    n = n(),
    ability_min = min(Ability),
    ability_max = max(Ability),
    delta_mean = mean(Delta_g),
    delta_median = median(Delta_g),
    delta_sd = sd(Delta_g),
    delta_se = delta_sd / sqrt(n),
    delta_ci_low_95 = delta_mean - qt(0.975, df = n - 1) * delta_se,
    delta_ci_high_95 = delta_mean + qt(0.975, df = n - 1) * delta_se,
    r2_mean = mean(Coherence_R2_1F),
    r2_sd = sd(Coherence_R2_1F),
    r2_se = r2_sd / sqrt(n),
    r2_ci_low_95 = r2_mean - qt(0.975, df = n - 1) * r2_se,
    r2_ci_high_95 = r2_mean + qt(0.975, df = n - 1) * r2_se,
    .groups = "drop"
  )

top_polluters <- metrics %>%
  arrange(desc(Delta_g)) %>%
  select(Model, Delta_g, Ability, Coherence_R2_1F, Coherence_Cosine_Raw, params_num, model_type) %>%
  slice_head(n = 30)

top_pillars <- metrics %>%
  arrange(Delta_g) %>%
  select(Model, Delta_g, Ability, Coherence_R2_1F, Coherence_Cosine_Raw, params_num, model_type) %>%
  slice_head(n = 30)

top_ability <- metrics %>%
  arrange(desc(Ability)) %>%
  select(Model, Ability, Delta_g, Coherence_R2_1F, params_num, model_type) %>%
  slice_head(n = 10)

bottom_ability <- metrics %>%
  arrange(Ability) %>%
  select(Model, Ability, Delta_g, Coherence_R2_1F, params_num, model_type) %>%
  slice_head(n = 10)

summary_tbl <- tibble(
  sample_models = N,
  sample_items = P,
  pc1_variance_pct = 100 * baseline_pc1_var,
  corr_ability_cosine_raw = cor(metrics$Ability, metrics$Coherence_Cosine_Raw),
  corr_ability_cosine_centered = cor(metrics$Ability, metrics$Coherence_Cosine_Centered),
  corr_ability_r2_1f = cor(metrics$Ability, metrics$Coherence_R2_1F),
  corr_ability_delta = cor(metrics$Ability, metrics$Delta_g),
  corr_delta_r2_1f = cor(metrics$Delta_g, metrics$Coherence_R2_1F),
  linear_r2 = summary(linear_fit)$r.squared,
  quad_r2 = summary(quad_fit)$r.squared,
  anova_linear_vs_quad_p = anova_lq$`Pr(>F)`[2],
  full_model_r2 = summary(full_fit)$r.squared,
  quad_turning_point = turning_point_est,
  quad_peak_delta = peak_delta_est
)

summary_tbl <- summary_tbl %>%
  bind_cols(
    boot_ci %>%
      transmute(
        corr_ability_delta_ci_low_95 = cor_ability_delta_lo,
        corr_ability_delta_ci_high_95 = cor_ability_delta_hi,
        corr_ability_r2_ci_low_95 = cor_ability_r2_lo,
        corr_ability_r2_ci_high_95 = cor_ability_r2_hi,
        corr_delta_r2_ci_low_95 = cor_delta_r2_lo,
        corr_delta_r2_ci_high_95 = cor_delta_r2_hi,
        linear_r2_ci_low_95 = linear_r2_lo,
        linear_r2_ci_high_95 = linear_r2_hi,
        quad_r2_ci_low_95 = quad_r2_lo,
        quad_r2_ci_high_95 = quad_r2_hi,
        full_r2_ci_low_95 = full_r2_lo,
        full_r2_ci_high_95 = full_r2_hi,
        quad_turning_point_ci_low_95 = turning_point_lo,
        quad_turning_point_ci_high_95 = turning_point_hi,
        quad_peak_delta_ci_low_95 = peak_delta_lo,
        quad_peak_delta_ci_high_95 = peak_delta_hi
      )
  )

all_regression_tbl <- bind_rows(
  linear_coef_tbl %>% mutate(model = "linear"),
  quad_coef_tbl %>% mutate(model = "quadratic"),
  full_coef_tbl %>% mutate(model = "full")
) %>%
  relocate(model, .before = term)

write_csv(metrics, file.path(outdir, "model_metrics.csv"))
write_csv(summary_tbl, file.path(outdir, "summary_core.csv"))
write_csv(decile_tbl, file.path(outdir, "ability_decile_summary.csv"))
write_csv(top_polluters, file.path(outdir, "top_polluters.csv"))
write_csv(top_pillars, file.path(outdir, "top_pillars.csv"))
write_csv(top_ability, file.path(outdir, "top_ability_models.csv"))
write_csv(bottom_ability, file.path(outdir, "bottom_ability_models.csv"))
write_csv(boot_tbl, file.path(outdir, "bootstrap_correlations.csv"))
write_csv(boot_ci, file.path(outdir, "bootstrap_correlation_ci.csv"))
write_csv(linear_coef_tbl, file.path(outdir, "delta_regression_linear_coefficients.csv"))
write_csv(quad_coef_tbl, file.path(outdir, "delta_regression_quadratic_coefficients.csv"))
write_csv(full_coef_tbl, file.path(outdir, "delta_regression_full_coefficients.csv"))
write_csv(all_regression_tbl, file.path(outdir, "delta_regression_all_coefficients.csv"))

png(file.path(outdir, "figure_scree_top20.png"), width = 1100, height = 700, res = 130)
plot(
  x = 1:20,
  y = 100 * var_pct[1:20],
  type = "b",
  pch = 19,
  xlab = "Principal Component",
  ylab = "Variance Explained (%)",
  main = "Scree Plot (Top 20 PCs)"
)
grid()
dev.off()

p1 <- ggplot(metrics, aes(x = Ability, y = Delta_g, color = model_type)) +
  geom_point(alpha = 0.65, size = 2.2) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "black", linewidth = 0.9) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  scale_color_manual(values = c(
    base_or_other = "#5c6b73",
    instruction_or_chat = "#4c78a8",
    merge_or_moe = "#f58518"
  )) +
  labs(
    x = "General Ability g (M = 100, SD = 15)",
    y = expression(Delta * "PC1 variance (jackknife impact)"),
    color = "Model Type",
    title = "Structural Impact vs Ability",
    subtitle = "Linear fit shown for correlation display; nonlinear modeling reported separately"
  ) +
  theme_minimal(base_size = 12)
ggsave(file.path(outdir, "figure_ability_vs_delta.png"), p1, width = 10, height = 6.5, dpi = 170)

p2 <- ggplot(metrics, aes(x = Ability, y = Coherence_R2_1F, color = model_type)) +
  geom_point(alpha = 0.6, size = 2.1) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "black", linewidth = 0.9) +
  scale_color_manual(values = c(
    base_or_other = "#5c6b73",
    instruction_or_chat = "#4c78a8",
    merge_or_moe = "#f58518"
  )) +
  labs(
    x = "General Ability g (M = 100, SD = 15)",
    y = expression("One-factor fit coherence (" * R^2 * ")"),
    color = "Model Type",
    title = "Structural Coherence vs Ability (Alternative Metric)"
  ) +
  theme_minimal(base_size = 12)
ggsave(file.path(outdir, "figure_ability_vs_r2.png"), p2, width = 10, height = 6.5, dpi = 170)

p3 <- ggplot(metrics, aes(x = Coherence_Cosine_Raw, y = Ability)) +
  geom_point(alpha = 0.6, size = 2, color = "#1f77b4") +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.9) +
  labs(
    x = "Raw cosine coherence",
    y = "General Ability g (M = 100, SD = 15)",
    title = "Diagnostic: Raw Cosine Coherence Is Nearly Collinear with Ability"
  ) +
  theme_minimal(base_size = 12)
ggsave(file.path(outdir, "figure_diagnostic_raw_cosine_collinearity.png"), p3, width = 9.5, height = 6, dpi = 170)

p4 <- ggplot(metrics, aes(x = Coherence_R2_1F, y = Delta_g, color = model_type)) +
  geom_point(alpha = 0.6, size = 2.1) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "black", linewidth = 0.9) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  scale_color_manual(values = c(
    base_or_other = "#5c6b73",
    instruction_or_chat = "#4c78a8",
    merge_or_moe = "#f58518"
  )) +
  labs(
    x = expression("One-factor fit coherence (" * R^2 * ")"),
    y = expression(Delta * "PC1 variance (jackknife impact)"),
    color = "Model Type",
    title = "Structural Impact vs Structural Coherence"
  ) +
  theme_minimal(base_size = 12)
ggsave(file.path(outdir, "figure_delta_vs_r2.png"), p4, width = 10, height = 6.5, dpi = 170)

critique_md <- c(
  "# Alignment-Tax Audit: Key Diagnostics",
  "",
  sprintf("- Sample after filtering: **%d models x %d items**.", N, P),
  sprintf("- PC1 variance explained: **%.2f%%**.", 100 * baseline_pc1_var),
  sprintf("- Raw cosine coherence vs ability correlation: **%.3f** (near-collinearity diagnostic).", cor(metrics$Ability, metrics$Coherence_Cosine_Raw)),
  sprintf("- Alternative coherence (1-factor R2) vs ability correlation: **%.3f**.", cor(metrics$Ability, metrics$Coherence_R2_1F)),
  sprintf("- Jackknife impact vs ability correlation: **%.3f**.", cor(metrics$Ability, metrics$Delta_g)),
  sprintf("- Jackknife impact vs 1-factor R2 correlation: **%.3f**.", cor(metrics$Delta_g, metrics$Coherence_R2_1F)),
  sprintf("- Quadratic fit significantly improves over linear for impact~ability: p = **%.3g**.", anova_lq$`Pr(>F)`[2]),
  "",
  "## Interpretation Guardrails",
  "",
  "- The original raw-cosine coherence metric is highly coupled to PC1 ability and should be treated as a diagnostic, not a standalone construct-validity endpoint.",
  "- The stronger supported pattern is a **middle-ability fracture** (inverted-U impact profile), not a monotonic 'smarter = more fractured' trend.",
  "- High-ability models in this snapshot are often near-neutral or negative impact; many strongest polluters cluster in the middle ability range.",
  "- Model-type heuristics (merge/instruction keywords) are noisy labels and should be interpreted as exploratory."
)
writeLines(critique_md, file.path(outdir, "analysis_critique.md"))

cat("Wrote outputs to:", outdir, "\n")
cat("Done.\n")
