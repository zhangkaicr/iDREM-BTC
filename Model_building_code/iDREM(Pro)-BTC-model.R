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

df2$age_grade <- factor(df2$age_grade, 
                         levels = c("younger", "older"),  # 第一个水平为参照组
                         ordered = FALSE)
df2$tumor_burden_grade <- factor(df2$tumor_burden_grade, 
                         levels = c("within", "beyond"),  # 第一个水平为参照组
                         ordered = FALSE)
df2$tumor_stage <- factor(df2$tumor_stage, 
                         levels = c("locally_advanced", "metastatic_intrahepatic", "metastatic_extrahepatic"),  # 第一个水平为参照组
                         ordered = FALSE)
df2$ecog_ps <- factor(df2$ecog_ps, 
                         levels = c("zero", "one", "two"),  # 第一个水平为参照组
                         ordered = FALSE)
df2$TBil_grade <- factor(df2$TBil_grade, 
                         levels = c("low", "high"),  # 第一个水平为参照组
                         ordered = FALSE)
df2$PD_L1 <- factor(df2$PD_L1, 
                         levels = c("low", "high"),  # 第一个水平为参照组
                         ordered = FALSE)
df2$MSI <- factor(df2$MSI, 
                         levels = c("low", "high"),  # 第一个水平为参照组
                         ordered = FALSE)
df2$TMB <- factor(df2$TMB, 
                         levels = c("low", "high"),  # 第一个水平为参照组
                         ordered = FALSE)


df$age_grade <- factor(df$age_grade, 
                         levels = c("younger", "older"),  # 第一个水平为参照组
                         ordered = FALSE)
df$tumor_burden_grade <- factor(df$tumor_burden_grade, 
                         levels = c("within", "beyond"),  # 第一个水平为参照组
                         ordered = FALSE)
df$tumor_stage <- factor(df$tumor_stage, 
                         levels = c("locally_advanced", "metastatic_intrahepatic", "metastatic_extrahepatic"),  # 第一个水平为参照组
                         ordered = FALSE)
df$ecog_ps <- factor(df$ecog_ps, 
                         levels = c("zero", "one", "two"),  # 第一个水平为参照组
                         ordered = FALSE)
df$TBil_grade <- factor(df$TBil_grade, 
                         levels = c("low", "high"),  # 第一个水平为参照组
                         ordered = FALSE)
df$PD_L1 <- factor(df$PD_L1, 
                         levels = c("low", "high"),  # 第一个水平为参照组
                         ordered = FALSE)
df$MSI <- factor(df$MSI, 
                         levels = c("low", "high"),  # 第一个水平为参照组
                         ordered = FALSE)
df$TMB <- factor(df$TMB, 
                         levels = c("low", "high"),  # 第一个水平为参照组
                         ordered = FALSE)


CoxFit2 <- coxph(Surv(months, status2) ~ age_grade + tumor_burden_grade + tumor_stage + ecog_ps + TBil_grade + CRP + CA199 + PD_L1 + MSI + TMB, data = df2)
tbl_regression(CoxFit2, exponentiate = TRUE)


#BEST
fm9 <- lme(LOG_CRP ~ ns(month, 3) + NLR, data = df,
           random = ~ ns(month, 3) | id, control = lmeControl(opt = 'optim'))
tab_model(fm9)


#BEST
fm21 <- lme(LOG_CA199 ~ ns(month, 3) + tumor_burden, data = df,
           random = ~ ns(month, 3) | id, control = lmeControl(opt = 'optim'))
tab_model(fm21)


#BEST
fm6 <- mixed_model(
    TBil_grade ~ month + tumor_stage,   
    data = df,
    random = ~ month | id,
    family = binomial(),  # 二项分布
)
tab_model(fm6)


jointFit2 <- 
  jm(Surv_object = CoxFit2, #生存cox方程
     Mixed_objects = list(fm9,fm21,fm6),  #纵向子模型/混合效应模型
     time_var = "month", #时间变量列名称
     id_var = "id", #患者id列名称
     n_thin = 1L, n_iter = 15000L, n_burnin = 5000L)
#查看结果
summary(jointFit2)


save.image("JMbayes2初步整体流程auc.RData")
load("JMbayes2初步整体流程auc.RData")
