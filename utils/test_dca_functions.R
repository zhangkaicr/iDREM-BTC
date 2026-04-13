# 设置脚本工作目录为当前文件所在目录，确保相对路径输出稳定
setwd("c:/Users/85330/Desktop/dynamic_dca")

# 输出当前R版本，便于和主脚本运行环境一致性核对
cat("R.version:", R.version.string, "\n")

# 统一结果目录：所有输出都写入 report 文件夹，不存在则自动创建
report_dir <- "c:/Users/85330/Desktop/dynamic_dca/report"
if (!dir.exists(report_dir)) dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

# 加载主函数脚本（会定义 dynamic_jm_dca 与 dynamic_jm_dca_batch）
source("c:/Users/85330/Desktop/dynamic_dca/dca.R")

# 定义一个统一的测试执行器：记录测试名、状态、错误信息和输出路径
run_case <- function(case_name, expr) {
  # 用 tryCatch 捕获错误，保证后续测试继续执行
  out <- tryCatch(
    {
      # 执行具体测试表达式
      res <- eval(expr)
      # 返回成功结构
      list(case = case_name, status = "ok", message = "", result = res)
    },
    error = function(e) {
      # 返回失败结构
      list(case = case_name, status = "failed", message = conditionMessage(e), result = NULL)
    }
  )
  # 打印当前测试结果摘要
  cat(sprintf("[CASE] %s -> %s\n", out$case, out$status))
  if (nzchar(out$message)) cat("  message:", out$message, "\n")
  # 返回结构化结果
  out
}

# 准备一个含二分类纵向变量的数据（ascites_bin: Yes=1, No=0）
long_bin <- dplyr::mutate(JMbayes2::pbc2, ascites_bin = as.integer(ascites == "Yes"))

# 统一的基础参数（减少重复代码）
base_args <- list(
  survival_data = JMbayes2::pbc2.id,
  survival_time_var = "years",
  survival_status_var = "status",
  visit_time_var = "year",
  cox_vars_str = "age + drug",
  id_var = "id",
  status_event_value = "dead",
  output_dir = report_dir
)

# 测试1：单个连续纵向子模型
case1 <- run_case(
  "single_longitudinal_model",
  quote(do.call(dynamic_jm_dca, c(
    list(
      horizon_time = 8,
      landmark_time = 5,
      longitudinal_data = JMbayes2::pbc2,
      lme_fixed = list(log(serBilir) ~ splines::ns(year, df = 2)),
      lme_random = list(~ splines::ns(year, df = 2) | id),
      output_prefix = "test_single"
    ),
    base_args
  )))
)

# 测试2：两个连续纵向子模型
case2 <- run_case(
  "two_longitudinal_models",
  quote(do.call(dynamic_jm_dca, c(
    list(
      horizon_time = 8,
      landmark_time = 5,
      longitudinal_data = JMbayes2::pbc2,
      lme_fixed = list(
        log(serBilir) ~ splines::ns(year, df = 2),
        albumin ~ year
      ),
      lme_random = list(
        ~ splines::ns(year, df = 2) | id,
        ~ 1 | id
      ),
      jm_control = list(n_chains = 1L, n_iter = 220L, n_burnin = 100L, n_thin = 2L),
      output_prefix = "test_two_models"
    ),
    base_args
  )))
)

# 测试3：三个连续纵向子模型
case3 <- run_case(
  "three_longitudinal_models",
  quote(do.call(dynamic_jm_dca, c(
    list(
      horizon_time = 8,
      landmark_time = 5,
      longitudinal_data = JMbayes2::pbc2,
      lme_fixed = list(
        log(serBilir) ~ splines::ns(year, df = 2),
        albumin ~ year,
        prothrombin ~ year
      ),
      lme_random = list(
        ~ splines::ns(year, df = 2) | id,
        ~ 1 | id,
        ~ 1 | id
      ),
      jm_control = list(n_chains = 1L, n_iter = 220L, n_burnin = 100L, n_thin = 2L),
      output_prefix = "test_three_models"
    ),
    base_args
  )))
)

# 测试4：二分类子模型（非严格检查，允许提示后映射）
case4 <- run_case(
  "binary_submodel_non_strict",
  quote(do.call(dynamic_jm_dca, c(
    list(
      horizon_time = 8,
      landmark_time = 5,
      longitudinal_data = long_bin,
      lme_fixed = list(
        log(serBilir) ~ splines::ns(year, df = 2),
        ascites_bin ~ year
      ),
      lme_random = list(
        ~ splines::ns(year, df = 2) | id,
        ~ 1 | id
      ),
      binary_submodel_index = c(2),
      binary_positive_value = 1,
      binary_check_strict = FALSE,
      jm_control = list(n_chains = 1L, n_iter = 220L, n_burnin = 100L, n_thin = 2L),
      output_prefix = "test_binary_non_strict"
    ),
    base_args
  )))
)

# 测试5：二分类严格检查（预期失败，用于确认检查逻辑生效）
case5 <- run_case(
  "binary_submodel_strict_expected_fail",
  quote(do.call(dynamic_jm_dca, c(
    list(
      horizon_time = 8,
      landmark_time = 5,
      longitudinal_data = JMbayes2::pbc2,
      lme_fixed = list(
        log(serBilir) ~ splines::ns(year, df = 2),
        ascites ~ year
      ),
      lme_random = list(
        ~ splines::ns(year, df = 2) | id,
        ~ 1 | id
      ),
      binary_submodel_index = c(2),
      binary_positive_value = "Yes",
      binary_check_strict = TRUE,
      output_prefix = "test_binary_strict"
    ),
    base_args
  )))
)

# 测试6：批量函数（pairwise，多组组合），体现如何批量分析
case6 <- run_case(
  "batch_pairwise_multi",
  quote(do.call(dynamic_jm_dca_batch, c(
    list(
      landmark_times = c(5, 6),
      horizon_times = c(8, 9),
      run_mode = "pairwise",
      continue_on_error = TRUE,
      longitudinal_data = JMbayes2::pbc2,
      lme_fixed = list(log(serBilir) ~ splines::ns(year, df = 2)),
      lme_random = list(~ splines::ns(year, df = 2) | id),
      jm_control = list(n_chains = 1L, n_iter = 220L, n_burnin = 100L, n_thin = 2L),
      output_prefix = "test_batch_pairwise"
    ),
    base_args
  )))
)

# 测试7：批量函数（grid 笛卡尔积），体现全组合批量分析
case7 <- run_case(
  "batch_grid_multi",
  quote(do.call(dynamic_jm_dca_batch, c(
    list(
      landmark_times = c(5, 5.5),
      horizon_times = c(8, 8.5),
      run_mode = "grid",
      continue_on_error = TRUE,
      longitudinal_data = JMbayes2::pbc2,
      lme_fixed = list(log(serBilir) ~ splines::ns(year, df = 2)),
      lme_random = list(~ splines::ns(year, df = 2) | id),
      jm_control = list(n_chains = 1L, n_iter = 220L, n_burnin = 100L, n_thin = 2L),
      output_prefix = "test_batch_grid"
    ),
    base_args
  )))
)

# 汇总所有测试结果
all_cases <- list(case1, case2, case3, case4, case5, case6, case7)
summary_df <- dplyr::bind_rows(lapply(all_cases, function(x) {
  data.frame(
    case = x$case,
    status = x$status,
    message = x$message,
    stringsAsFactors = FALSE
  )
}))

# 将汇总结果写出为CSV，便于后续留档
readr::write_csv(summary_df, file.path(report_dir, "test_summary_dca_functions.csv"))

# 打印测试汇总
cat("\n===== TEST SUMMARY =====\n")
print(summary_df)

# 额外展示 batch 函数返回的 summary，明确批量分析输出结构
if (!is.null(case6$result)) {
  cat("\n===== BATCH PAIRWISE SUMMARY =====\n")
  print(case6$result$summary)
  cat("pairwise n_ok:", case6$result$n_ok, " n_failed:", case6$result$n_failed, "\n")
}
if (!is.null(case7$result)) {
  cat("\n===== BATCH GRID SUMMARY =====\n")
  print(case7$result$summary)
  cat("grid n_ok:", case7$result$n_ok, " n_failed:", case7$result$n_failed, "\n")
}

# 如果有非预期失败（除 strict_expected_fail 外），返回非零状态便于CI识别
unexpected_failed <- summary_df$status == "failed" & summary_df$case != "binary_submodel_strict_expected_fail"
if (any(unexpected_failed)) {
  quit(status = 1)
}
