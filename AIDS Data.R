# =========================================================
# 最终完整且可运行的脚本
# R Script for GEE and WGEE Analysis of AIDS Data
# =========================================================

# --- 1. 加载库与数据 (Load Libraries & Data) ---
# 确保在运行前已经安装了所有必要的包
# install.packages("geepack")
# install.packages("survival")
# install.packages("survminer")
# install.packages("dplyr")
# install.packages("ggplot2")
# install.packages("tibble")
# install.packages("forcats")
#
# 请确保wgeesel已经从本地文件安装
# install.packages("C:/Users/Lenovo/Desktop/wgeesel_1.5.tar.gz", repos = NULL, type = "source")

library(wgeesel)
library(ggplot2)
library(JM)
library(geepack)
library(dplyr)
library(tibble)
library(forcats)
library(survival)
library(survminer)


# --- 2. 加载并重构数据 (Load and Reshape Data) ---
data(aids)

# 将数据从不平衡的长格式转换为平衡格式
N <- length(unique(aids$patient))
m <- 5 # 固定的时间点数量
Y=id=drug1=gender1=prevOI1=AZT1=matrix(NA, N, m)

for(i in 1:N){
  id[i,] <- rep(i, m)
  # 填充已有的CD4观测值
  Y[i, 1:length(aids$CD4[aids$patient==i])] <- aids$CD4[aids$patient==i]
  # 填充基线协变量
  drug1[i,] <- rep(aids$drug[aids$patient==i][1], m)
  gender1[i,] <- rep(aids$gender[aids$patient==i][1], m)
  prevOI1[i,] <- rep(aids$prevOI[aids$patient==i][1], m)
  AZT1[i,] <- rep(aids$AZT[aids$patient==i][1], m)
}
Timee <- rep(c(0, 2, 6, 12, 18), N)

# 创建最终的分析数据集
Data <- data.frame(id=as.numeric(t(id)), Time=Timee, CD4=as.numeric(t(Y)),
                   drug=as.numeric(t(drug1))-1,
                   gender=as.numeric(t(gender1))-1, prevOI=as.numeric(t(prevOI1))-1,
                   AZT=as.numeric(t(AZT1))-1)
# 创建缺失指示变量 R
R <- rep(1, N*m)
R[is.na(Data$CD4)==TRUE] <- 0
Data$R <- R


# --- 3. 拟合两个模型 (Fit Models) ---

# WGEE 模型 (处理MAR缺失)
fit <- wgee(CD4 ~ Time + Time*drug + gender + prevOI + AZT, 
            data=Data, id=Data$id, family="gaussian",
            corstr="exchangeable", scale=NULL, 
            mismodel = R ~ Time + drug + prevOI + gender)
summary(fit)

# 标准 GEE 模型 (仅使用完整观测)
Data_complete <- na.omit(Data)
fit_gee <- geeglm(
  CD4 ~ Time + Time*drug + gender + prevOI + AZT,
  data    = Data_complete,
  id      = Data_complete$id, # 【修正点】: 明确使用Data_complete中的id向量
  family  = gaussian(),
  corstr  = "exchangeable"
)
summary(fit_gee)


# --- 4. 准备画图所需的数据 (Prepare Data for Plotting) ---

# 定义模型公式 (方便复用)
gee_formula <- CD4 ~ Time + Time*drug + gender + prevOI + AZT

# 计算 GEE 拟合值
Data_complete$mu_gee <- as.numeric(
  predict(fit_gee, newdata = Data_complete, type = "response")
)

# 计算 WGEE 拟合值 (手动进行矩阵乘法 X * beta)
Xcc <- model.matrix(gee_formula, data = Data_complete)
# 【修正点】: 直接将系数矩阵转为向量，依靠自然顺序进行计算
b_wgee_vec <- as.vector(fit$beta) 
Data_complete$mu_wgee <- as.numeric(Xcc %*% b_wgee_vec)

# 计算残差
Data_complete$res_gee  <- Data_complete$CD4 - Data_complete$mu_gee
Data_complete$res_wgee <- Data_complete$CD4 - Data_complete$mu_wgee


# --- 5. 开始绘制所有图表 (Generate All Plots) ---

# (图0) 探索性分析：个体CD4轨迹图 (Spaghetti Plot)
p0 <- ggplot(aids, aes(x = obstime, y = CD4, group = patient)) +
    geom_line(alpha = 0.5, aes(color = drug)) +  # 按药物给线条上色
    facet_wrap(~ drug) + # 按药物分面
    scale_x_continuous(breaks = c(0, 2, 6, 12, 18)) +
    scale_color_manual(values = c("#1b9e77", "#d95f02")) + 
    labs(title = "Longitudinal CD4 Cell Count Trajectories by Treatment Group",
         x = "Observation Time (weeks)",
         y = "CD4 Count") +
    theme_minimal() +
    theme(legend.position = "none")


ggsave(
  filename = "C:/Users/Lenovo/Desktop/plot_0_CD4_Trajectories.pdf",
  plot = p0,    # 你之前定义的图对象
  width = 10,   # 宽度（英寸）
  height = 8,   # 高度（英寸）
  units = "in"  # 单位
)







# (图1) 生存曲线 (Kaplan-Meier Survival Curve)
surv_obj <- Surv(time = aids$Time, event = aids$death)
fit_km <- survfit(surv_obj ~ drug, data = aids)
library(cowplot)

# ggsurvplot 对象
p_km <- ggsurvplot(
  fit_km, data = aids,
  conf.int = TRUE, pval = TRUE, risk.table = TRUE,
  ggtheme = theme_minimal(),
  palette = c("#1b9e77", "#d95f02")
)

# cowplot 拼接主图和风险表
km_combined <- plot_grid(
  p_km$plot,      # 主图
  p_km$table,     # 风险表
  ncol = 1,
  rel_heights = c(3, 1)  # 调整主图和风险表比例
)

# 保存 PDF
ggsave(
  filename = "C:/Users/Lenovo/Desktop/plot_KM_withRisk.pdf",
  plot = km_combined,
  width = 10, height = 8, units = "in"
)








# (图2) 基线CD4分布 (Baseline CD4 Distribution)
baseline_df <- Data %>%
  filter(Time == 0) %>%
  mutate(drug = factor(drug, levels = c(0,1), labels = c("ddC","ddI")))
# (图2) 基线CD4分布
p2 <- ggplot(baseline_df, aes(x = CD4)) +
  geom_histogram(aes(y = after_stat(density)), bins = 20, alpha = 0.6, fill="grey80", color = "white") +
  geom_density(aes(color = drug), linewidth = 1.2) +
  facet_wrap(~ drug, ncol = 2, scales = "free_y") +
  scale_color_manual(values = c("ddC" = "#1b9e77", "ddI" = "#d95f02")) +
  labs(title = "Baseline CD4 Distribution by Drug Group",
       x = "Baseline CD4 Count", y = "Density") +
  theme_minimal() +
  theme(legend.position = "none")

# 保存到桌面
ggsave(
  filename = "C:/Users/Lenovo/Desktop/plot_2_BaselineCD4.pdf",
  plot = p2,
  width = 8, height = 6, units = "in"
)






# (图3) 观测 vs 预测均值图 (Observed vs. Predicted Mean Trajectories)
plot_df <- Data_complete %>%
  mutate(drug = factor(drug, levels = c(0,1), labels = c("ddC","ddI"))) %>%
  group_by(Time, drug) %>%
  summarise(
    obs_mean  = mean(CD4, na.rm = TRUE),
    obs_se    = sd(CD4, na.rm = TRUE) / sqrt(n()),
    pred_wgee = mean(mu_wgee, na.rm = TRUE), # 添加na.rm以防万一
    pred_gee  = mean(mu_gee, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(lo = obs_mean - 1.96*obs_se, hi = obs_mean + 1.96*obs_se)

# (图3) 观测 vs 预测均值图 (Observed vs. Predicted Mean Trajectories)
p3 <- ggplot(plot_df, aes(x = Time)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, fill = "grey70") +
  geom_point(aes(y = obs_mean), color = "black") +
  geom_line(aes(y = pred_wgee, color = "WGEE"), linetype = "dashed", linewidth=1) +
  geom_line(aes(y = pred_gee,  color = "GEE"), linetype = "solid", linewidth=1) +
  facet_wrap(~ drug) +
  scale_color_manual(name = "Model", values = c("GEE" = "#1f77b4", "WGEE" = "#d62728")) +
  labs(title = "Observed vs Predicted Mean CD4 by Drug", 
       x = "Weeks", y = "CD4 Count") +
  theme_minimal() + 
  theme(legend.position = "top")

# 保存到桌面
ggsave(
  filename = "C:/Users/Lenovo/Desktop/plot_3_Obs_vs_Pred.pdf",
  plot = p3,
  width = 8, height = 6, units = "in"
)






# (图4) 系数对比图 (Coefficient Comparison Plot)
sum_gee <- summary(fit_gee)
sum_wgee <- summary(fit)
gee_results <- sum_gee$coefficients %>% as.data.frame() %>% rownames_to_column(var = "variable") %>%
  transmute(model = "GEE", variable, estimate = Estimate, std_error = `Std.err`,
            conf.low = estimate - 1.96 * std_error, conf.high = estimate + 1.96 * std_error)
wgee_results <- data.frame(model="WGEE", variable=sum_wgee$coefnames, estimate=as.numeric(sum_wgee$beta),
                           std_error=as.numeric(sum_wgee$se.robust)) %>%
  mutate(conf.low=estimate - 1.96*std_error, conf.high=estimate + 1.96*std_error)
coef_df <- bind_rows(gee_results, wgee_results) %>%
  mutate(model = factor(model, levels = c("WGEE", "GEE")))

# (图4) 系数对比图 (Coefficient Comparison Plot)
p4 <- ggplot(coef_df, aes(x = estimate, y = fct_rev(variable), color = model)) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), 
                 height = 0.2, position = position_dodge(width = 0.5), linewidth = 0.8) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = c("GEE" = "#1f77b4", "WGEE" = "#d62728")) +
  labs(
    title = "Comparison of GEE and WGEE Model Coefficients", 
    subtitle = "Points are estimates; Bars are 95% confidence intervals",
    x = "Coefficient Estimate", 
    y = "Predictor Variable", 
    color = "Model"
  ) +
  theme_minimal() + 
  theme(legend.position = "top")

# 保存到桌面
ggsave(
  filename = "C:/Users/Lenovo/Desktop/plot_4_Coefficient_Comparison.pdf",
  plot = p4,
  width = 9, height = 6, units = "in"
)














# (图5) 失访概率图 (Dropout Probability Plot)
mis_fit <- fit$mis_fit
pred_data <- expand.grid(Time = seq(0, 18, length.out = 50), drug = c(0, 1),
                         prevOI = c(0, 1), gender = 0)
pred_data$prob_stay <- predict(mis_fit, newdata = pred_data, type = "response")
pred_data <- pred_data %>%
  mutate(drug = factor(drug, levels=c(0,1), labels=c("ddC","ddI")),
         prevOI = factor(prevOI, levels=c(0,1), labels=c("No Previous OI","Previous OI")))

# (图5) 失访概率图 (Dropout Probability Plot)
p5 <- ggplot(pred_data, aes(x=Time, y=prob_stay, color=prevOI, linetype=prevOI)) +
  geom_line(linewidth = 1.2) + 
  facet_wrap(~ drug) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  labs(
    title = "Probability of Remaining in Study (Not Dropping Out)", 
    subtitle = "Based on the WGEE missingness model",
    x = "Observation Time (weeks)", 
    y = "Probability of Observation", 
    color = "Patient History", 
    linetype = "Patient History"
  ) +
  theme_minimal() + 
  theme(legend.position = "top")

# 保存到桌面
ggsave(
  filename = "C:/Users/Lenovo/Desktop/plot_5_Dropout_Probability.pdf",
  plot = p5,
  width = 9, height = 6, units = "in"
)










# --- 6. 模型诊断图 (Model Diagnostic Plots) ---
# 创建一个用于诊断图的长格式数据框
res_long <- bind_rows(
  transmute(Data_complete, model = "GEE",  fitted = mu_gee,  resid = res_gee),
  transmute(Data_complete, model = "WGEE", fitted = mu_wgee, resid = res_wgee)
) %>%
  filter(is.finite(fitted), is.finite(resid)) %>%
  mutate(model = factor(model, levels = c("GEE","WGEE")))

# (图6) 残差 vs 拟合值图 (Residuals vs. Fitted Plot)
p6 <- ggplot(res_long, aes(fitted, resid)) +
  geom_point(alpha = 0.45, size = 0.9) +
  geom_smooth(method = "loess", se = FALSE, color = "red") +
  facet_wrap(~ model, nrow = 1) +
  labs(title = "Residuals vs Fitted (Left: GEE, Right: WGEE)",
       x = "Fitted values", y = "Residuals") +
  theme_minimal()

ggsave(
  filename = "C:/Users/Lenovo/Desktop/plot_6_Residuals_vs_Fitted.pdf",
  plot = p6,
  width = 9, height = 6, units = "in"
)


# (图7) 正态 Q-Q 图 (Normal Q-Q Plot)
p7 <- ggplot(res_long, aes(sample = resid)) +
  stat_qq() + 
  stat_qq_line() +
  facet_wrap(~ model, nrow = 1) +
  labs(title = "Normal Q-Q Plot of Residuals (Left: GEE, Right: WGEE)",
       x = "Theoretical Quantiles", y = "Sample Quantiles") +
  theme_minimal()

ggsave(
  filename = "C:/Users/Lenovo/Desktop/plot_7_QQ.pdf",
  plot = p7,
  width = 9, height = 6, units = "in"
)
