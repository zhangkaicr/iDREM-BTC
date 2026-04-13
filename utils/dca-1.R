# 加载并静默启动所需 R 包，避免控制台输出冗余信息
suppressPackageStartupMessages({
  library(dplyr)      # 数据清洗与变换
  library(stringr)    # 字符串匹配与处理
  library(readr)      # 高效读写 CSV
  library(ggplot2)    # 可视化
  library(JMbayes2) # 联合模型拟合与预测
  library(nlme)     # 线性混合效应模型
  library(survival)   # Cox 生存模型
  library(dcurves)  # DCA（决策曲线分析）
})

rm(list = ls())
# 设置随机种子，保证结果可复现
set.seed(2026)

# 定义 Landmark 与 Horizon 时间点（单位：年）
landmark_time <- 5      # 预测起点：患者存活至 5 年时才进入分析
horizon_time <- 8       # 预测终点：预测 8 年时是否发生事件
threshold_grid <- seq(0.01, 0.50, by = 0.01)  # DCA 阈值网格（1%–50%）

# 读取 JMbayes2 自带 PBC 数据：纵向重复测量与基线生存信息
data_long <- as_tibble(JMbayes2::pbc2)   # 纵向数据（多次随访记录）
data_surv <- as_tibble(JMbayes2::pbc2.id) %>% 
  mutate(
    event = as.integer(status == "dead"),
    event_horizon = as.integer(event == 1L & years > landmark_time & years <= horizon_time)
  )

# 筛选在 landmark_time 仍存活的患者 ID，用于后续动态预测
ids_at_risk <- data_surv %>%
  filter(years > landmark_time) %>%
  pull(id) %>% unique()

# 构造 Landmark 数据集：仅保留上述患者 5 年及以前的测量记录
newdata_landmark <- data_long %>%
  filter(id %in% ids_at_risk, year <= landmark_time) %>%
  mutate(
    years = year,   # 与生存模型时间变量同名，便于后续合并
    event = 0L      # Landmark 时未发生事件
  )

# 拟合纵向混合效应模型：log(血清胆红素) 随时间变化采用 2 自由度样条
lme_fit <- nlme::lme(
  fixed = log(serBilir) ~ splines::ns(year, df = 2),
  random = ~ splines::ns(year, df = 2) | id,  # 允许个体随机截距与斜率
  data = as.data.frame(data_long),
  na.action = na.omit,
  control = nlme::lmeControl(opt = "optim")
)

# 拟合 Cox 生存模型：以年龄与药物分组为协变量
cox_fit <- survival::coxph(
  formula = survival::Surv(years, event) ~ age + drug,
  data = as.data.frame(data_surv),
  x = TRUE  # 保存设计矩阵，供 JMbayes2 调用
)

# 联合建模：将纵向 lme 与生存 cox 结合，评估生物标志物动态对生存的关联
jm_fit <- JMbayes2::jm(
  Surv_object = cox_fit,
  Mixed_objects = lme_fit,
  time_var = "year",
  n_chains = 1L,   # 链数（示例用 1 链，快速演示）
  n_iter = 500L,   # 总迭代
  n_burnin = 200L, # 丢弃前 200 次
  n_thin = 2L      # 每 2 次取 1 样本，降低自相关
)

# 对 Landmark 人群进行动态预测：得到 8 年生存概率
pred_surv_raw <- predict(
  object = jm_fit,
  newdata = as.data.frame(newdata_landmark),
  process = "event",      # 预测事件（生存）概率
  times = horizon_time,   # 预测时间点
  return_newdata = TRUE,  # 将预测结果合并回 newdata
  control = list(
    n_samples = 100L, # MC 样本数
    n_mcmc = 20L,     # 每条链额外 MCMC 步
    parallel = "snow",
    cores = 1L
  )
)

# 将预测结果转为 tibble，并自动识别生存概率列
# 将 predict() 返回的对象强制转换为 tibble，方便后续管道操作
pred_surv_tbl <- as_tibble(pred_surv_raw)

# 记录原始 landmark 数据框的所有列名，用于后续排除
source_cols <- names(newdata_landmark)

# 从预测结果中剔除原始列，仅保留 predict() 新生成的列（即预测相关列）
candidate_cols <- setdiff(names(pred_surv_tbl), source_cols)

# 在候选列中进一步筛选：仅保留列名包含 pred/surv/event/risk 等关键词的列（忽略大小写）
candidate_cols <- candidate_cols[str_detect(candidate_cols, regex("pred|surv|event|risk", ignore_case = TRUE))]

# 在上一步候选列中，再排除带 low/upp/ci/sd/se 等字样（通常代表置信区间或标准误）的列，优先保留点估计列
priority_cols <- candidate_cols[!str_detect(candidate_cols, regex("low|upp|ci|sd|se", ignore_case = TRUE))]

# 优先使用 priority_cols 的第一个有效列；若为空，则退而求其次使用 candidate_cols 的第一个有效列
pred_col <- dplyr::coalesce(first(priority_cols), first(candidate_cols))

# 如果最终仍无法识别到合适的预测列，则抛出错误并打印所有可用列名，方便调试
if (is.na(pred_col)) {
  stop(sprintf("未识别到JMbayes2事件预测列，请检查predict返回字段: %s", paste(names(pred_surv_tbl), collapse = ", ")))
}

# 构造 DCA 所需数据框：每人保留最后一次测量，计算事件概率
# 从预测结果 tibble 中构造 DCA 分析所需的最简数据框：每人仅保留最后一次随访记录，并计算事件概率
dca_input <- pred_surv_tbl %>%                       # 以预测结果 tibble 为起点
  select(id, year, !!sym(pred_col)) %>%               # 仅保留患者 ID、随访时间及预测生存概率列（pred_col 为动态识别的列名）
  group_by(id) %>%                                    # 按患者分组，后续每组只取一条记录
  slice_max(order_by = year, n = 1, with_ties = FALSE) %>%  # 取每组最晚一次随访记录（Landmark 最近记录）
  ungroup() %>%                                       # 解除分组，恢复普通 tibble
  rename(pred_raw = !!sym(pred_col)) %>%
  left_join(select(data_surv, id, event_horizon), by = "id") %>%  # 把生存结局（8 年内是否事件）拼接到预测结果上
  mutate(pred_event = if (str_detect(pred_col, regex("surv", ignore_case = TRUE))) 1 - pred_raw else pred_raw) %>%
  filter(is.finite(pred_event), pred_event >= 0, pred_event <= 1) %>%  # 剔除非法概率值（NA、负值或>1）
  select(id, event_horizon, pred_event)               # 最终仅保留三列：患者 ID、实际事件标志、预测事件概率

# 执行决策曲线分析：评估不同阈值下的净收益
dca_fit <- dcurves::dca(
  formula = event_horizon ~ pred_event,
  data = as.data.frame(dca_input),
  thresholds = threshold_grid,
  label = list(pred_event = "JMbayes2动态预测")
)

# 绘制 DCA 曲线：对比“全部干预”与“全不干预”策略
dca_plot <- plot(dca_fit, smooth = TRUE) +
  labs(
    title = "基于JMbayes2动态预测的DCA曲线",
    subtitle = sprintf("Landmark = %.1f 年, Horizon = %.1f 年", landmark_time, horizon_time)
  )

# 保存图像与数据：便于后续报告与复现
ggsave(
  filename = "dca_curve_jmbayes2.png",
  plot = dca_plot,
  width = 8,
  height = 6,
  dpi = 300
)

write_csv(dca_input, "dca_input_from_dynamic_prediction.csv")
saveRDS(dca_fit, "dca_fit_jmbayes2.rds")
saveRDS(jm_fit, "jm_fit_jmbayes2.rds")
