library(JMbayes2) # 用于构建动态预测模型，包含联合建模和动态预测功能
library(survival) # 用于生存分析，提供生存曲线、Cox回归等功能
library(tidyverse) # 数据科学工具集，包含数据处理、可视化等功能
library(dlookr) # 用于探索性数据分析，提供数据诊断和统计分析功能
library(visdat) # 用于探索性数据分析的可视化，可快速查看数据结构和缺失值
library(DT) # 用于创建交互式表格，使数据展示更加整洁美观
library(gtsummary) # 用于生成标准化的基线特征表，支持统计检验
library(sjPlot) # 用于统计可视化和数据展示，提供多种图表类型
library(ggpubr) # 基于ggplot2的统计图形包，提供发表级别的图形功能
# 项目管理包,用于管理项目文件路径
library(here)
# 生存分析包,用于生存分析相关函数
# Excel文件读写包
library(openxlsx)
# 数据处理和可视化包集合,包含dplyr,ggplot2等
library(tidyverse)
# 用于制作统计表格和模型输出的包
library(finalfit)
library(survminer)

# 2.导入数据-------------------
rm(list = ls())  # 清空工作环境，避免变量冲突

# 读取数据并查看数据结构
list.files("input")  # 查看数据目录下的文件
df <- read.xlsx(here("input","df_long.xlsx"))  # 读取Excel文件中的cofounders工作表
glimpse(df)  # 查看数据框的结构，包括变量类型和前几行数据


# 读取数据并查看数据结构
list.files("input")  # 查看数据目录下的文件
df2 <- read.xlsx(here("input","df.xlsx"))  # 读取Excel文件中的cofounders工作表
glimpse(df2)  # 查看数据框的结构，包括变量类型和前几行数据


# 1. 首先确保df_male数据的预处理
df <- df %>%
    mutate(
        # 因子型变量
        age_grade = as.factor(age_grade),
        tumor_burden_grade = as.factor(tumor_burden_grade),
        tumor_stage = as.factor(tumor_stage),
        ecog_ps = as.factor(ecog_ps),
        TBil_grade = as.factor(TBil_grade),
        id = as.factor(id),
        
        # 数值型变量
        CRP = as.numeric(CRP),
        CA199 = as.numeric(CA199),
        LOG_CRP = as.numeric(LOG_CRP),
        NLR = as.numeric(NLR),
        LOG_CA199 = as.numeric(LOG_CA199),
        tumor_burden = as.numeric(tumor_burden),
        month = as.numeric(month),
        months = as.numeric(months),
        status2 = as.numeric(status2)
    )

set.seed(123)
test_data <- slice_sample(.data = df2,
                         prop = 0.3, 
                         by = "status2")

# 从pbc2.id数据集中筛选出不在test_data中的样本作为训练集
# 使用setdiff()函数找出pbc2.id中不在test_data中的id
# filter()函数基于id进行筛选,保留这些id对应的样本
train_data <- 
  df2 |> 
  filter(id %in% setdiff(df2$id,test_data$id))


train_data_long <- 
  df |> 
  filter(id %in% train_data$id) |> 
  mutate(team = "train")

test_data_long <- 
  df |> 
  filter(id %in% test_data$id) |> 
  mutate(team = "test")


data_all_long <- 
  train_data_long |> 
  bind_rows(test_data_long)




train_data$age_grade <- factor(train_data$age_grade, 
                         levels = c("younger", "older"),  # 第一个水平为参照组
                         ordered = FALSE)
train_data$tumor_burden_grade <- factor(train_data$tumor_burden_grade, 
                         levels = c("within", "beyond"),  # 第一个水平为参照组
                         ordered = FALSE)
train_data$tumor_stage <- factor(train_data$tumor_stage, 
                         levels = c("locally_advanced", "metastatic_intrahepatic", "metastatic_extrahepatic"),  # 第一个水平为参照组
                         ordered = FALSE)
train_data$ecog_ps <- factor(train_data$ecog_ps, 
                         levels = c("zero", "one", "two"),  # 第一个水平为参照组
                         ordered = FALSE)
train_data$TBil_grade <- factor(train_data$TBil_grade, 
                         levels = c("low", "high"),  # 第一个水平为参照组
                         ordered = FALSE)



train_data_long$age_grade <- factor(train_data_long$age_grade, 
                         levels = c("younger", "older"),  # 第一个水平为参照组
                         ordered = FALSE)
train_data_long$tumor_burden_grade <- factor(train_data_long$tumor_burden_grade, 
                         levels = c("within", "beyond"),  # 第一个水平为参照组
                         ordered = FALSE)
train_data_long$tumor_stage <- factor(train_data_long$tumor_stage, 
                         levels = c("locally_advanced", "metastatic_intrahepatic", "metastatic_extrahepatic"),  # 第一个水平为参照组
                         ordered = FALSE)
train_data_long$ecog_ps <- factor(train_data_long$ecog_ps, 
                         levels = c("zero", "one", "two"),  # 第一个水平为参照组
                         ordered = FALSE)
train_data_long$TBil_grade <- factor(train_data_long$TBil_grade, 
                         levels = c("low", "high"),  # 第一个水平为参照组
                         ordered = FALSE)






CoxFit1 <- coxph(Surv(months, status2) ~ age_grade + tumor_burden_grade + tumor_stage + ecog_ps + TBil_grade + CRP + CA199, data = train_data)
tbl_regression(CoxFit1, exponentiate = TRUE)



#BEST
fm9 <- lme(LOG_CRP ~ ns(month, 3) + NLR, data = train_data_long,
           random = ~ ns(month, 3) | id, control = lmeControl(opt = 'optim'))
tab_model(fm9)


#BEST
fm21 <- lme(LOG_CA199 ~ ns(month, 3) + tumor_burden, data = train_data_long,
           random = ~ ns(month, 3) | id, control = lmeControl(opt = 'optim'))
tab_model(fm21)


#BEST
fm6 <- mixed_model(
    TBil_grade ~ month + tumor_stage,   
    data = train_data_long,
    random = ~ month | id,
    family = binomial(),  # 二项分布
)
tab_model(fm6)




jointFit1 <- 
  jm(Surv_object = CoxFit1, #生存cox方程
     Mixed_objects = list(fm9,fm21,fm6),  #纵向子模型/混合效应模型
     time_var = "month", #时间变量列名称
     id_var = "id", #患者id列名称
     n_thin = 1L, n_iter = 12000L, n_burnin = 2000L)
#查看结果
summary(jointFit1)

save.image("JMbayes2初步整体流程subgroup.RData")
load("JMbayes2初步整体流程subgroup.RData")
