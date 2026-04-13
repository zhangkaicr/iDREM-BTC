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
  binary_submodel_index = integer(0),
  binary_positive_value = 1L,
  binary_check_strict = FALSE,
  binary_family = stats::binomial(),
  id_var = "id",
  status_event_value = 1L,
  lme_control = list(opt = "optim"),
  glmm_control = list(),
  jm_control = list(n_chains = 1L, n_iter = 500L, n_burnin = 200L, n_thin = 2L),
  predict_control = list(n_samples = 100L, n_mcmc = 20L, parallel = "snow", cores = 1L),
  threshold_grid = seq(0.01, 0.50, by = 0.01),
  output_dir = "report",
  output_prefix = "dca_curve",
  smooth = TRUE
) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
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

  lme_fixed_list <- if (inherits(lme_fixed, "formula")) list(lme_fixed) else lme_fixed
  lme_random_list <- if (inherits(lme_random, "formula")) list(lme_random) else lme_random
  if (!is.list(lme_fixed_list) || !is.list(lme_random_list)) {
    stop("lme_fixed 和 lme_random 需为公式或公式列表。")
  }
  if (length(lme_fixed_list) != length(lme_random_list)) {
    stop("lme_fixed 与 lme_random 的长度必须一致。")
  }
  if (length(lme_fixed_list) < 1 || length(lme_fixed_list) > 3) {
    stop("当前函数仅支持 1 到 3 个纵向子模型。")
  }
  binary_submodel_index <- sort(unique(as.integer(binary_submodel_index)))
  if (any(is.na(binary_submodel_index)) || any(binary_submodel_index < 1) || any(binary_submodel_index > length(lme_fixed_list))) {
    stop("binary_submodel_index 必须是 1 到纵向子模型数量之间的整数索引。")
  }

  binary_positive_map <- setNames(vector("list", length(binary_submodel_index)), as.character(binary_submodel_index))
  if (length(binary_submodel_index) > 0) {
    if (is.list(binary_positive_value)) {
      if (length(binary_positive_value) != length(binary_submodel_index)) {
        stop("当 binary_positive_value 为列表时，其长度必须与 binary_submodel_index 一致。")
      }
      binary_positive_map <- setNames(binary_positive_value, as.character(binary_submodel_index))
    } else {
      binary_positive_map <- setNames(rep(list(binary_positive_value), length(binary_submodel_index)), as.character(binary_submodel_index))
    }
  }

  long_vars_required <- unique(c(
    id_var,
    visit_time_var,
    unlist(lapply(lme_fixed_list, all.vars)),
    unlist(lapply(lme_random_list, all.vars))
  ))
  missing_long_vars <- setdiff(long_vars_required, names(long_tbl))
  if (length(missing_long_vars) > 0) {
    stop(sprintf("纵向数据缺少变量: %s", paste(missing_long_vars, collapse = ", ")))
  }
  long_model_data <- long_tbl %>%
    filter(if_all(all_of(long_vars_required), ~ !is.na(.x)))

  if (length(binary_submodel_index) > 0) {
    for (idx in binary_submodel_index) {
      response_var <- all.vars(lme_fixed_list[[idx]])[1]
      positive_value <- binary_positive_map[[as.character(idx)]]
      response_vec <- long_model_data[[response_var]]
      raw_unique <- sort(unique(response_vec[!is.na(response_vec)]))
      if (binary_check_strict) {
        if (!(length(raw_unique) == 2 && identical(raw_unique, c(0, 1)))) {
          stop(sprintf(
            "二分类子模型索引 %s 的结果变量 '%s' 不是严格0/1二分类，检测到取值: %s",
            idx,
            response_var,
            paste(raw_unique, collapse = ", ")
          ))
        }
      } else if (!(length(raw_unique) == 2 && identical(raw_unique, c(0, 1)))) {
        message(sprintf(
          "提示: 二分类子模型索引 %s 的结果变量 '%s' 原始取值不是严格0/1，检测到: %s；将按 binary_positive_value 映射。",
          idx,
          response_var,
          paste(raw_unique, collapse = ", ")
        ))
      }
      long_model_data[[response_var]] <- if (is.logical(response_vec)) {
        as.integer(response_vec)
      } else if (is.numeric(response_vec) || is.integer(response_vec)) {
        as.integer(response_vec %in% positive_value)
      } else {
        as.integer(as.character(response_vec) %in% as.character(positive_value))
      }
    }
  }

  newdata_landmark <- long_model_data %>%
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

  lme_control_list <- if (!is.list(lme_control) || is.null(names(lme_control))) {
    rep(list(lme_control), length(lme_fixed_list))
  } else if (all(vapply(lme_control, is.list, logical(1)))) {
    if (length(lme_control) != length(lme_fixed_list)) {
      stop("当 lme_control 提供列表列表时，其长度必须与纵向子模型数量一致。")
    }
    lme_control
  } else {
    rep(list(lme_control), length(lme_fixed_list))
  }

  mixed_objects <- Map(
    f = function(fixed_formula, random_formula, control_list_item, model_index) {
      fit_data <- as.data.frame(long_model_data)
      if (model_index %in% binary_submodel_index) {
        GLMMadaptive::mixed_model(
          fixed = fixed_formula,
          random = random_formula,
          data = fit_data,
          family = binary_family,
          na.action = na.omit,
          control = glmm_control
        )
      } else {
        nlme::lme(
          fixed = fixed_formula,
          random = random_formula,
          data = fit_data,
          na.action = na.omit,
          control = do.call(nlme::lmeControl, control_list_item)
        )
      }
    },
    fixed_formula = lme_fixed_list,
    random_formula = lme_random_list,
    control_list_item = lme_control_list,
    model_index = seq_along(lme_fixed_list)
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
        Mixed_objects = if (length(mixed_objects) == 1) mixed_objects[[1]] else mixed_objects,
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
    mixed_objects = mixed_objects,
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

dynamic_jm_dca_batch <- function(
  landmark_times,
  horizon_times,
  run_mode = c("grid", "pairwise"),
  continue_on_error = TRUE,
  fail_if_all_failed = TRUE,
  ...
) {
  run_mode <- match.arg(run_mode)
  if (run_mode == "pairwise") {
    if (length(landmark_times) != length(horizon_times)) {
      stop("run_mode='pairwise' 时 landmark_times 与 horizon_times 长度必须一致。")
    }
    combos <- data.frame(landmark_time = landmark_times, horizon_time = horizon_times)
  } else {
    combos <- expand.grid(landmark_time = landmark_times, horizon_time = horizon_times)
  }

  # 预分配结果容器，保持与组合行号一一对应
  results <- vector("list", nrow(combos))
  summary_tbl <- vector("list", nrow(combos))

  for (i in seq_len(nrow(combos))) {
    lm_t <- combos$landmark_time[i]
    hz_t <- combos$horizon_time[i]
    # 捕获每一组参数的运行结果，失败时返回 error 对象
    run_res <- tryCatch(
      dynamic_jm_dca(landmark_time = lm_t, horizon_time = hz_t, ...),
      error = function(e) e
    )
    if (inherits(run_res, "error")) {
      if (!continue_on_error) stop(run_res$message)
      # 用 [i] <- list(NULL) 保留索引位，避免 [[i]] <- NULL 删除元素
      results[i] <- list(NULL)
      summary_tbl[[i]] <- data.frame(
        landmark_time = lm_t,
        horizon_time = hz_t,
        status = "failed",
        message = run_res$message,
        png = NA_character_,
        pdf = NA_character_
      )
    } else {
      # 成功结果按索引位写入，保证 results 与 summary 行号一致
      results[i] <- list(run_res)
      summary_tbl[[i]] <- data.frame(
        landmark_time = lm_t,
        horizon_time = hz_t,
        status = "ok",
        message = "",
        png = run_res$paths$png,
        pdf = run_res$paths$pdf
      )
    }
  }

  summary_df <- dplyr::bind_rows(summary_tbl)
  # 统计批量运行整体成功/失败数量，便于上游脚本快速判断
  n_failed <- sum(summary_df$status == "failed")
  n_ok <- sum(summary_df$status == "ok")
  if (fail_if_all_failed && n_ok == 0) {
    stop("批量运行未产生任何成功结果；请检查 summary$message 并适当提高 jm_control 迭代参数。")
  }
  list(
    summary = summary_df,
    results = results,
    n_ok = n_ok,
    n_failed = n_failed
  )
}
