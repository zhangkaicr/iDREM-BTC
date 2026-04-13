suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(readr)
  library(ggplot2)
  library(JMbayes2)
  library(nlme)
  library(survival)
  library(dcurves)
})

set.seed(2026)

dynamic_jm_dca <- function(
  horizon_time,
  landmark_time,
  longitudinal_data,
  survival_data,
  survival_time_var,
  survival_status_var,
  visit_time_var,
  cox_vars_str,
  lme_fixed,
  lme_random,
  id_var = "id",
  status_event_value = 1L,
  lme_control = list(opt = "optim"),
  jm_control = list(n_chains = 1L, n_iter = 500L, n_burnin = 200L, n_thin = 2L),
  predict_control = list(n_samples = 100L, n_mcmc = 20L, parallel = "snow", cores = 1L),
  threshold_grid = seq(0.01, 0.50, by = 0.01),
  output_dir = ".",
  output_prefix = "dca_curve",
  smooth = TRUE
) {
  long_tbl <- as_tibble(longitudinal_data)
  surv_tbl <- as_tibble(survival_data)

  status_vec <- surv_tbl[[survival_status_var]]
  event_binary <- if (is.logical(status_vec)) {
    as.integer(status_vec)
  } else if (is.numeric(status_vec) || is.integer(status_vec)) {
    as.integer(status_vec %in% status_event_value)
  } else {
    as.integer(as.character(status_vec) %in% as.character(status_event_value))
  }

  surv_tbl <- surv_tbl %>%
    mutate(
      .surv_time = .data[[survival_time_var]],
      .event_binary = event_binary,
      event_horizon = as.integer(.event_binary == 1L & .surv_time > landmark_time & .surv_time <= horizon_time)
    )

  ids_at_risk <- surv_tbl %>%
    filter(.surv_time > landmark_time) %>%
    pull(all_of(id_var)) %>%
    unique()

  cox_covars <- all.vars(as.formula(paste("~", cox_vars_str)))

  newdata_landmark <- long_tbl %>%
    filter(.data[[id_var]] %in% ids_at_risk, .data[[visit_time_var]] <= landmark_time) %>%
    mutate(
      .surv_time = .data[[visit_time_var]],
      .event_binary = 0L
    )

  missing_covars <- setdiff(cox_covars, names(newdata_landmark))
  if (length(missing_covars) > 0) {
    newdata_landmark <- newdata_landmark %>%
      left_join(
        surv_tbl %>%
          select(all_of(c(id_var, missing_covars))),
        by = id_var
      )
  }

  lme_fit <- nlme::lme(
    fixed = lme_fixed,
    random = lme_random,
    data = as.data.frame(long_tbl),
    na.action = na.omit,
    control = do.call(nlme::lmeControl, lme_control)
  )

  cox_formula <- as.formula(sprintf("survival::Surv(.surv_time, .event_binary) ~ %s", cox_vars_str))
  cox_fit <- survival::coxph(
    formula = cox_formula,
    data = as.data.frame(surv_tbl),
    x = TRUE
  )

  jm_fit <- do.call(
    JMbayes2::jm,
    c(
      list(
        Surv_object = cox_fit,
        Mixed_objects = lme_fit,
        time_var = visit_time_var
      ),
      jm_control
    )
  )

  pred_surv_raw <- predict(
    object = jm_fit,
    newdata = as.data.frame(newdata_landmark),
    process = "event",
    times = horizon_time,
    return_newdata = TRUE,
    control = predict_control
  )

  pred_surv_tbl <- as_tibble(pred_surv_raw)
  source_cols <- names(newdata_landmark)
  candidate_cols <- setdiff(names(pred_surv_tbl), source_cols)
  candidate_cols <- candidate_cols[str_detect(candidate_cols, regex("pred|surv|event|risk", ignore_case = TRUE))]
  priority_cols <- candidate_cols[!str_detect(candidate_cols, regex("low|upp|ci|sd|se", ignore_case = TRUE))]
  pred_col <- dplyr::coalesce(first(priority_cols), first(candidate_cols))
  if (is.na(pred_col)) {
    stop(sprintf("未识别到JMbayes2事件预测列，请检查predict返回字段: %s", paste(names(pred_surv_tbl), collapse = ", ")))
  }

  cox_lp_tbl <- surv_tbl %>%
    filter(.data[[id_var]] %in% ids_at_risk) %>%
    mutate(
      cox_lp = as.numeric(predict(cox_fit, newdata = as.data.frame(.), type = "lp")),
      cox_lp_prob = plogis(cox_lp)
    ) %>%
    select(all_of(id_var), cox_lp, cox_lp_prob)

  dca_input <- pred_surv_tbl %>%
    select(all_of(id_var), all_of(visit_time_var), !!sym(pred_col)) %>%
    group_by(.data[[id_var]]) %>%
    slice_max(order_by = .data[[visit_time_var]], n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    rename(pred_raw = !!sym(pred_col)) %>%
    left_join(
      surv_tbl %>% select(all_of(id_var), event_horizon),
      by = id_var
    ) %>%
    mutate(pred_event = if (str_detect(pred_col, regex("surv", ignore_case = TRUE))) 1 - pred_raw else pred_raw) %>%
    left_join(cox_lp_tbl, by = id_var) %>%
    filter(
      is.finite(pred_event), pred_event >= 0, pred_event <= 1,
      is.finite(cox_lp_prob), cox_lp_prob >= 0, cox_lp_prob <= 1
    ) %>%
    select(all_of(id_var), event_horizon, pred_event, cox_lp, cox_lp_prob)

  dca_fit <- dcurves::dca(
    formula = event_horizon ~ pred_event + cox_lp_prob,
    data = as.data.frame(dca_input),
    thresholds = threshold_grid,
    as_probability = c("pred_event", "cox_lp_prob"),
    label = list(
      pred_event = "JMbayes2动态死亡风险",
      cox_lp_prob = "Cox线性预测值映射风险"
    )
  )

  dca_plot <- plot(dca_fit, smooth = smooth) +
    labs(
      title = "JMbayes2与Cox线性预测值的DCA比较",
      subtitle = sprintf("Landmark = %.2f, Horizon = %.2f", landmark_time, horizon_time)
    )

  lm_tag <- gsub("\\.", "p", format(landmark_time, trim = TRUE))
  hz_tag <- gsub("\\.", "p", format(horizon_time, trim = TRUE))
  file_stem <- sprintf("%s_LM%s_HZ%s", output_prefix, lm_tag, hz_tag)
  png_path <- file.path(output_dir, paste0(file_stem, ".png"))
  pdf_path <- file.path(output_dir, paste0(file_stem, ".pdf"))
  csv_path <- file.path(output_dir, paste0(file_stem, "_input.csv"))
  dca_rds_path <- file.path(output_dir, paste0(file_stem, "_dca.rds"))
  jm_rds_path <- file.path(output_dir, paste0(file_stem, "_jm.rds"))

  ggsave(filename = png_path, plot = dca_plot, width = 8, height = 6, dpi = 300)
  ggsave(filename = pdf_path, plot = dca_plot, width = 8, height = 6)
  readr::write_csv(dca_input, csv_path)
  saveRDS(dca_fit, dca_rds_path)
  saveRDS(jm_fit, jm_rds_path)

  list(
    dca_input = dca_input,
    dca_fit = dca_fit,
    dca_plot = dca_plot,
    jm_fit = jm_fit,
    cox_fit = cox_fit,
    paths = list(
      png = png_path,
      pdf = pdf_path,
      input_csv = csv_path,
      dca_rds = dca_rds_path,
      jm_rds = jm_rds_path
    )
  )
}

result_example <- dynamic_jm_dca(
  horizon_time = 8,
  landmark_time = 5,
  longitudinal_data = JMbayes2::pbc2,
  survival_data = JMbayes2::pbc2.id,
  survival_time_var = "years",
  survival_status_var = "status",
  visit_time_var = "year",
  cox_vars_str = "age + drug",
  lme_fixed = log(serBilir) ~ splines::ns(year, df = 2),
  lme_random = ~ splines::ns(year, df = 2) | id,
  id_var = "id",
  status_event_value = "dead",
  output_dir = ".",
  output_prefix = "dca_curve_jmbayes2"
)

print(result_example$paths)
