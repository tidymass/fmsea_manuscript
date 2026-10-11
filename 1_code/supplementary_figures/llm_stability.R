library(r4projects)
setwd(get_project_wd())
rm(list = ls())

library(featuremsea)
library(tidyverse)

out_dir <- "3_data_analysis/9_llm_stability"
raw_dir <- file.path(out_dir, "raw")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

n_runs <- 10
sample_sources <- c("plasma", "urine")

# Fig. 2c input and the scores shown in Fig. 2c (single GPT-5.5 run per label)
load("2_data/7_pregnancy_plasma/plasma_res.rda")
load("2_data/7_pregnancy_plasma/urine_res.rda")

reported_scores <- bind_rows(
  plasma_res@significant_modules %>%
    transmute(sample_source = "plasma", pathway_id, reported_score = matrix_confidence_score),
  urine_res@significant_modules %>%
    transmute(sample_source = "urine", pathway_id, reported_score = matrix_confidence_score)
)

models <- list(
  `GPT-5.5` = list(provider = "openai", model = "gpt-5.5",
                   api_key = Sys.getenv("OPENAI_API_KEY")),
  `Qwen3-32B` = list(provider = "siliconflow", model = "Qwen/Qwen3-32B",
                     api_key = Sys.getenv("SILICONFLOW_API_KEY")),
  `Qwen3.6-35B-A3B` = list(provider = "siliconflow", model = "Qwen/Qwen3.6-35B-A3B",
                           api_key = Sys.getenv("SILICONFLOW_API_KEY"))
)

# analyze_matrix_relevance() only needs an S4 object with a significant_modules slot
setClass("llm_input", representation(significant_modules = "data.frame"))

llm_input <- new(
  "llm_input",
  significant_modules = plasma_res@significant_modules %>%
    dplyr::select(pathway_id, pathway_name, pathway_description) %>%
    as.data.frame()
)
rm(plasma_res, urine_res)

#####-------------------------------------------------------------------------
##### run LLM scoring (each run is cached, so the loop can be resumed)
#####-------------------------------------------------------------------------
run_once <- function(model_name, sample_source, run, max_attempt = 3) {
  file <- file.path(raw_dir, paste0("Fig2c_", model_name, "_", sample_source, "_run", run, ".rds"))
  if (file.exists(file)) {
    return(readRDS(file))
  }

  m <- models[[model_name]]
  for (attempt in seq_len(max_attempt)) {
    res <- tryCatch(
      analyze_matrix_relevance(
        results = llm_input,
        sample_source = sample_source,
        api_key = m$api_key,
        provider = m$provider,
        model = m$model
      )@significant_modules,
      error = function(e) {
        message("  failed: ", conditionMessage(e))
        NULL
      }
    )
    if (!is.null(res) && any(!is.na(res$matrix_confidence_score))) break
    res <- NULL
  }
  if (is.null(res)) {
    warning(model_name, " / ", sample_source, " / run ", run, " failed after ", max_attempt, " attempts")
    return(NULL)
  }

  res <- res %>%
    transmute(model = model_name, sample_source = sample_source, run = run,
              pathway_id, pathway_name, score = matrix_confidence_score,
              reason = matrix_confidence_reason)
  saveRDS(res, file)
  res
}

run_grid <- expand_grid(model = names(models), sample_source = sample_sources, run = seq_len(n_runs))

llm_scores <- pmap_dfr(run_grid, function(model, sample_source, run) {
  message(model, " | ", sample_source, " | run ", run)
  run_once(model, sample_source, run)
})

save(llm_scores, file = file.path(out_dir, "llm_scores.rda"), compress = "xz")

#####-------------------------------------------------------------------------
##### agreement statistics
#####-------------------------------------------------------------------------
anchors <- c(0, 25, 50, 75, 100)

# most frequent score; ties are broken by closeness to the median
mode_score <- function(x) {
  x <- x[!is.na(x)]
  tab <- table(x)
  cand <- as.numeric(names(tab)[tab == max(tab)])
  cand[which.min(abs(cand - median(x)))]
}

# Krippendorff's alpha; rows = units (metabolite sets), columns = raters (runs)
kripp_alpha <- function(mat, level = c("ordinal", "interval")) {
  level <- match.arg(level)
  mat <- mat[rowSums(!is.na(mat)) >= 2, , drop = FALSE]
  vals <- sort(unique(as.vector(mat[!is.na(mat)])))
  if (length(vals) < 2) return(NA_real_)
  k <- length(vals)
  o <- matrix(0, k, k)
  for (u in seq_len(nrow(mat))) {
    x <- mat[u, !is.na(mat[u, ])]
    n_uc <- tabulate(match(x, vals), nbins = k)
    o <- o + (outer(n_uc, n_uc) - diag(n_uc)) / (length(x) - 1)
  }
  n_c <- rowSums(o)
  n <- sum(n_c)
  if (level == "interval") {
    d2 <- outer(vals, vals, function(a, b) (a - b)^2)
  } else {
    d2 <- matrix(0, k, k)
    for (i in seq_len(k)) for (j in seq_len(k)) {
      lo <- min(i, j)
      hi <- max(i, j)
      d2[i, j] <- (sum(n_c[lo:hi]) - (n_c[i] + n_c[j]) / 2)^2
    }
  }
  1 - (n - 1) * sum(o * d2) / sum(outer(n_c, n_c) * d2)
}

# quadratic weighted Cohen's kappa on the fixed 5-level scale
weighted_kappa <- function(x, y) {
  keep <- !is.na(x) & !is.na(y)
  x <- factor(x[keep], levels = anchors)
  y <- factor(y[keep], levels = anchors)
  obs <- table(x, y) / sum(keep)
  expd <- outer(rowSums(obs), colSums(obs))
  w <- outer(seq_along(anchors), seq_along(anchors), function(i, j) (i - j)^2)
  1 - sum(w * obs) / sum(w * expd)
}

paired_wilcox_p <- function(x, y) {
  tryCatch(
    suppressWarnings(wilcox.test(x, y, paired = TRUE, exact = FALSE)$p.value),
    error = function(e) NA_real_
  )
}

# 1. within-model stability per metabolite set and matrix label
set_stability <- llm_scores %>%
  group_by(model, sample_source, pathway_id, pathway_name) %>%
  summarise(
    n_valid = sum(!is.na(score)),
    mode_score = mode_score(score),
    median_score = median(score, na.rm = TRUE),
    mode_agreement = mean(score == mode_score, na.rm = TRUE),
    n_distinct_scores = n_distinct(score[!is.na(score)]),
    score_range = diff(range(score, na.rm = TRUE)),
    score_sd = sd(score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(reported_scores, by = c("sample_source", "pathway_id"))

# rows = metabolite set (x matrix label, if not a grouping variable), columns = runs
score_matrix <- function(df) {
  df %>%
    dplyr::select(any_of("sample_source"), pathway_id, run, score) %>%
    pivot_wider(names_from = run, values_from = score) %>%
    dplyr::select(-any_of("sample_source"), -pathway_id) %>%
    as.matrix()
}

summarise_stability <- function(d, s) {
  tibble(
    n_sets = nrow(s),
    n_runs = n_distinct(d$run),
    parse_failure_rate = mean(is.na(d$score)),
    mean_mode_agreement = mean(s$mode_agreement),
    prop_sets_identical = mean(s$n_distinct_scores == 1),
    prop_sets_within_one_anchor = mean(s$score_range <= 25),
    max_score_range = max(s$score_range),
    kripp_alpha_ordinal = kripp_alpha(score_matrix(d), "ordinal"),
    kripp_alpha_interval = kripp_alpha(score_matrix(d), "interval")
  )
}

within_model_summary <- bind_rows(
  llm_scores %>%
    group_by(model, sample_source) %>%
    group_modify(function(d, g) {
      summarise_stability(d, set_stability %>% semi_join(g, by = c("model", "sample_source")))
    }) %>%
    ungroup(),
  llm_scores %>%
    group_by(model) %>%
    group_modify(function(d, g) {
      summarise_stability(d, set_stability %>% semi_join(g, by = "model"))
    }) %>%
    ungroup() %>%
    mutate(sample_source = "both")
) %>%
  arrange(model, sample_source)

# agreement of the scores shown in Fig. 2c with the GPT-5.5 modal scores
reported_vs_mode <- set_stability %>%
  filter(model == "GPT-5.5") %>%
  group_by(sample_source) %>%
  summarise(
    n_sets = n(),
    prop_reported_equal_mode = mean(reported_score == mode_score, na.rm = TRUE),
    prop_reported_within_one_anchor = mean(abs(reported_score - mode_score) <= 25, na.rm = TRUE),
    .groups = "drop"
  )

# 2. cross-model agreement on consensus (modal) scores, for each pair of models
consensus <- set_stability %>%
  dplyr::select(sample_source, pathway_id, pathway_name, model, mode_score) %>%
  pivot_wider(names_from = model, values_from = mode_score)

summarise_cross <- function(x, y) {
  tibble(
    n_items = length(x),
    exact_agreement = mean(x == y, na.rm = TRUE),
    within_one_anchor = mean(abs(x - y) <= 25, na.rm = TRUE),
    weighted_kappa = weighted_kappa(x, y),
    spearman_rho = suppressWarnings(cor(x, y, method = "spearman", use = "complete.obs")),
    kripp_alpha_ordinal = kripp_alpha(cbind(x, y), "ordinal"),
    mean_diff_model1_minus_model2 = mean(x - y, na.rm = TRUE),
    wilcoxon_p = paired_wilcox_p(x, y)
  )
}

cross_model_summary <- combn(names(models), 2, simplify = FALSE) %>%
  map_dfr(function(p) {
    bind_rows(
      consensus %>%
        group_by(sample_source) %>%
        group_modify(~ summarise_cross(.x[[p[1]]], .x[[p[2]]])) %>%
        ungroup(),
      summarise_cross(consensus[[p[1]]], consensus[[p[2]]]) %>%
        mutate(sample_source = "both")
    ) %>%
      mutate(model1 = p[1], model2 = p[2], .before = 1)
  })

# run-level pairwise weighted kappa over the 18 set x label items,
# within each model and between each pair of models
run_wide <- llm_scores %>%
  mutate(rater = paste(model, run, sep = "_")) %>%
  dplyr::select(sample_source, pathway_id, rater, score) %>%
  pivot_wider(names_from = rater, values_from = score)

raters <- setdiff(names(run_wide), c("sample_source", "pathway_id"))
pairwise_kappa <- combn(raters, 2, simplify = FALSE) %>%
  map_dfr(function(p) {
    m1 <- sub("_\\d+$", "", p[1])
    m2 <- sub("_\\d+$", "", p[2])
    tibble(
      rater1 = p[1], rater2 = p[2],
      comparison = if (m1 == m2) paste("Within", m1) else paste(m1, "vs", m2),
      weighted_kappa = weighted_kappa(run_wide[[p[1]]], run_wide[[p[2]]]),
      exact_agreement = mean(run_wide[[p[1]]] == run_wide[[p[2]]], na.rm = TRUE)
    )
  })

pairwise_kappa_summary <- pairwise_kappa %>%
  group_by(comparison) %>%
  summarise(
    n_pairs = n(),
    median_kappa = median(weighted_kappa),
    min_kappa = min(weighted_kappa),
    max_kappa = max(weighted_kappa),
    median_exact_agreement = median(exact_agreement),
    .groups = "drop"
  )

# 3. reproducibility of the Fig. 2c conclusion (plasma scores > urine scores)
# plasma and urine calls are independent; run r of plasma is paired with run r of urine
label_diff <- llm_scores %>%
  dplyr::select(model, run, sample_source, pathway_id, pathway_name, score) %>%
  pivot_wider(names_from = sample_source, values_from = score) %>%
  mutate(diff = plasma - urine)

label_diff_by_run <- label_diff %>%
  group_by(model, run) %>%
  summarise(
    n_plasma_higher = sum(diff > 0, na.rm = TRUE),
    n_equal = sum(diff == 0, na.rm = TRUE),
    n_urine_higher = sum(diff < 0, na.rm = TRUE),
    mean_diff = mean(diff, na.rm = TRUE),
    wilcoxon_p = paired_wilcox_p(plasma, urine),
    .groups = "drop"
  )

label_diff_by_set <- label_diff %>%
  group_by(model, pathway_id, pathway_name) %>%
  summarise(
    prop_plasma_higher = mean(diff > 0, na.rm = TRUE),
    prop_equal = mean(diff == 0, na.rm = TRUE),
    prop_urine_higher = mean(diff < 0, na.rm = TRUE),
    median_diff = median(diff, na.rm = TRUE),
    .groups = "drop"
  )

# the same comparison on modal scores, and on the scores shown in Fig. 2c
modal_diff <- set_stability %>%
  dplyr::select(model, sample_source, pathway_id, mode_score) %>%
  pivot_wider(names_from = sample_source, values_from = mode_score)

reported_diff <- reported_scores %>%
  pivot_wider(names_from = sample_source, values_from = reported_score) %>%
  mutate(model = "GPT-5.5 (Fig. 2c)")

summarise_label_diff <- function(d) {
  tibble(
    n_sets = nrow(d),
    n_plasma_higher = sum(d$plasma > d$urine, na.rm = TRUE),
    n_equal = sum(d$plasma == d$urine, na.rm = TRUE),
    n_urine_higher = sum(d$plasma < d$urine, na.rm = TRUE),
    mean_diff = mean(d$plasma - d$urine, na.rm = TRUE),
    wilcoxon_p = paired_wilcox_p(d$plasma, d$urine)
  )
}

# modal scores of each model vs the single run shown in Fig. 2c
label_diff_summary <- bind_rows(modal_diff, reported_diff) %>%
  group_by(model) %>%
  group_modify(~ summarise_label_diff(.x)) %>%
  ungroup()

# how often the Fig. 2c conclusion holds across the 10 runs of each model
label_diff_run_summary <- label_diff_by_run %>%
  group_by(model) %>%
  summarise(
    n_runs = n(),
    prop_runs_plasma_higher_majority = mean(n_plasma_higher > n_urine_higher),
    prop_runs_no_urine_higher = mean(n_urine_higher == 0),
    median_n_plasma_higher = median(n_plasma_higher),
    prop_runs_wilcoxon_p_lt_0.05 = mean(wilcoxon_p < 0.05, na.rm = TRUE),
    .groups = "drop"
  )

save(llm_scores, set_stability, within_model_summary, reported_vs_mode, consensus,
     cross_model_summary, pairwise_kappa, pairwise_kappa_summary,
     label_diff, label_diff_by_run, label_diff_by_set, label_diff_summary,
     label_diff_run_summary,
     file = file.path(out_dir, "llm_stability_results.rda"), compress = "xz")

write_csv(set_stability, file.path(out_dir, "set_stability.csv"))
write_csv(within_model_summary, file.path(out_dir, "within_model_summary.csv"))
write_csv(reported_vs_mode, file.path(out_dir, "reported_vs_mode.csv"))
write_csv(consensus, file.path(out_dir, "cross_model_consensus.csv"))
write_csv(cross_model_summary, file.path(out_dir, "cross_model_summary.csv"))
write_csv(pairwise_kappa_summary, file.path(out_dir, "pairwise_kappa_summary.csv"))
write_csv(label_diff_by_run, file.path(out_dir, "label_diff_by_run.csv"))
write_csv(label_diff_by_set, file.path(out_dir, "label_diff_by_set.csv"))
write_csv(label_diff_summary, file.path(out_dir, "label_diff_summary.csv"))
write_csv(label_diff_run_summary, file.path(out_dir, "label_diff_run_summary.csv"))

print(within_model_summary, width = Inf)
print(reported_vs_mode, width = Inf)
print(cross_model_summary, width = Inf)
print(pairwise_kappa_summary, width = Inf)
print(label_diff_summary, width = Inf)
print(label_diff_run_summary, width = Inf)
