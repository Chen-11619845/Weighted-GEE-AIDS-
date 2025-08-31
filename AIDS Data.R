# =========================================================
# --- 1. Load Libraries & Data ---
# Ensure all necessary packages are installed before running
# install.packages("geepack")
# install.packages("survival")
# install.packages("survminer")
# install.packages("dplyr")
# install.packages("ggplot2")
# install.packages("tibble")
# install.packages("forcats")


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


# --- 2. Load and Reshape Data ---
data(aids)

# Convert data from an unbalanced long format to a balanced format
N <- length(unique(aids$patient))
m <- 5 # Fixed number of time points
Y=id=drug1=gender1=prevOI1=AZT1=matrix(NA, N, m)

for(i in 1:N){
  id[i,] <- rep(i, m)
  # Fill in existing CD4 observations
  Y[i, 1:length(aids$CD4[aids$patient==i])] <- aids$CD4[aids$patient==i]
  # Fill in baseline covariates
  drug1[i,] <- rep(aids$drug[aids$patient==i][1], m)
  gender1[i,] <- rep(aids$gender[aids$patient==i][1], m)
  prevOI1[i,] <- rep(aids$prevOI[aids$patient==i][1], m)
  AZT1[i,] <- rep(aids$AZT[aids$patient==i][1], m)
}
Timee <- rep(c(0, 2, 6, 12, 18), N)

# Create the final analysis dataset
Data <- data.frame(id=as.numeric(t(id)), Time=Timee, CD4=as.numeric(t(Y)),
                   drug=as.numeric(t(drug1))-1,
                   gender=as.numeric(t(gender1))-1, prevOI=as.numeric(t(prevOI1))-1,
                   AZT=as.numeric(t(AZT1))-1)
# Create the missingness indicator variable R
R <- rep(1, N*m)
R[is.na(Data$CD4)==TRUE] <- 0
Data$R <- R


# --- 3. Fit Models ---

# WGEE Model (handles MAR missingness)
fit <- wgee(CD4 ~ Time + Time*drug + gender + prevOI + AZT, 
            data=Data, id=Data$id, family="gaussian",
            corstr="exchangeable", scale=NULL, 
            mismodel = R ~ Time + drug + prevOI + gender)
summary(fit)

# Standard GEE Model (uses only complete observations)
Data_complete <- na.omit(Data)
fit_gee <- geeglm(
  CD4 ~ Time + Time*drug + gender + prevOI + AZT,
  data    = Data_complete,
  id      = Data_complete$id, # [Correction]: Explicitly use the id vector from Data_complete
  family  = gaussian(),
  corstr  = "exchangeable"
)
summary(fit_gee)


# --- 4. Prepare Data for Plotting ---

# Define model formula (for reusability)
gee_formula <- CD4 ~ Time + Time*drug + gender + prevOI + AZT

# Calculate GEE fitted values
Data_complete$mu_gee <- as.numeric(
  predict(fit_gee, newdata = Data_complete, type = "response")
)

# Calculate WGEE fitted values (manual matrix multiplication X * beta)
Xcc <- model.matrix(gee_formula, data = Data_complete)
# [Correction]: Directly convert the coefficient matrix to a vector for calculation based on natural order
b_wgee_vec <- as.vector(fit$beta) 
Data_complete$mu_wgee <- as.numeric(Xcc %*% b_wgee_vec)

# Calculate residuals
Data_complete$res_gee  <- Data_complete$CD4 - Data_complete$mu_gee
Data_complete$res_wgee <- Data_complete$CD4 - Data_complete$mu_wgee


# --- 5. Generate All Plots ---

# (Plot 0) Exploratory Analysis: Individual CD4 Trajectories (Spaghetti Plot)
p0 <- ggplot(aids, aes(x = obstime, y = CD4, group = patient)) +
    geom_line(alpha = 0.5, aes(color = drug)) +  # Color lines by drug
    facet_wrap(~ drug) + # Facet by drug
    scale_x_continuous(breaks = c(0, 2, 6, 12, 18)) +
    scale_color_manual(values = c("#1b9e77", "#d95f02")) + 
    labs(title = "Longitudinal CD4 Cell Count Trajectories by Treatment Group",
         x = "Observation Time (weeks)",
         y = "CD4 Count") +
    theme_minimal() +
    theme(legend.position = "none")


ggsave(
  filename = "C:/Users/Lenovo/Desktop/plot_0_CD4_Trajectories.pdf",
  plot = p0,    # The plot object you defined earlier
  width = 10,   # width (inches)
  height = 8,   # height (inches)
  units = "in"  # units
)



# ggsurvplot object
p_km <- ggsurvplot(
  fit_km, data = aids,
  conf.int = TRUE, pval = TRUE, risk.table = TRUE,
  ggtheme = theme_minimal(),
  palette = c("#1b9e77", "#d95f02")
)

# Combine the main plot and risk table using cowplot
km_combined <- plot_grid(
  p_km$plot,      # main plot
  p_km$table,     # risk table
  ncol = 1,
  rel_heights = c(3, 1)  # Adjust the ratio of the main plot and risk table
)

# Save as PDF
ggsave(
  filename = "C:/Users/Lenovo/Desktop/plot_KM_withRisk.pdf",
  plot = km_combined,
  width = 10, height = 8, units = "in"
)








# (Plot 2) Baseline CD4 Distribution
baseline_df <- Data %>%
  filter(Time == 0) %>%
  mutate(drug = factor(drug, levels = c(0,1), labels = c("ddC","ddI")))

p2 <- ggplot(baseline_df, aes(x = CD4)) +
  geom_histogram(aes(y = after_stat(density)), bins = 20, alpha = 0.6, fill="grey80", color = "white") +
  geom_density(aes(color = drug), linewidth = 1.2) +
  facet_wrap(~ drug, ncol = 2, scales = "free_y") +
  scale_color_manual(values = c("ddC" = "#1b9e77", "ddI" = "#d95f02")) +
  labs(title = "Baseline CD4 Distribution by Drug Group",
       x = "Baseline CD4 Count", y = "Density") +
  theme_minimal() +
  theme(legend.position = "none")

# Save to Desktop
ggsave(
  filename = "C:/Users/Lenovo/Desktop/plot_2_BaselineCD4.pdf",
  plot = p2,
  width = 8, height = 6, units = "in"
)






# (Plot 3) Observed vs. Predicted Mean Trajectories
plot_df <- Data_complete %>%
  mutate(drug = factor(drug, levels = c(0,1), labels = c("ddC","ddI"))) %>%
  group_by(Time, drug) %>%
  summarise(
    obs_mean  = mean(CD4, na.rm = TRUE),
    obs_se    = sd(CD4, na.rm = TRUE) / sqrt(n()),
    pred_wgee = mean(mu_wgee, na.rm = TRUE), # Add na.rm just in case
    pred_gee  = mean(mu_gee, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(lo = obs_mean - 1.96*obs_se, hi = obs_mean + 1.96*obs_se)


p3 <- ggplot(plot_df, aes(x = Time)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, fill = "grey70") +
  geom_point(aes(y = obs_mean), color = "black") +
  geom_line(aes(y = pred_wgee, color = "WGEE"), linetype = "dashed", linewidth=1) +
  geom_line(aes(y = pred_gee,  color = "GEE"), linetype = "solid", linewidth=1) +
  facet_wrap(~ drug) +
  scale_color_manual(name = "Model", values = c("GEE" = "#1f77b4", "WGEE" = "#d62728")) +
  labs(title = "Observed vs Predicted Mean CD4 by Drug", 
       x = "Weeks", y = "CD4 Count") +
  theme_minimal() + 
  theme(legend.position = "top")

# Save to Desktop
ggsave(
  filename = "C:/Users/Lenovo/Desktop/plot_3_Obs_vs_Pred.pdf",
  plot = p3,
  width = 8, height = 6, units = "in"
)






# (Plot 4) Coefficient Comparison Plot
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

# Save to Desktop
ggsave(
  filename = "C:/Users/Lenovo/Desktop/plot_4_Coefficient_Comparison.pdf",
  plot = p4,
  width = 9, height = 6, units = "in"
)














# (Plot 5) Dropout Probability Plot
mis_fit <- fit$mis_fit
pred_data <- expand.grid(Time = seq(0, 18, length.out = 50), drug = c(0, 1),
                         prevOI = c(0, 1), gender = 0)
pred_data$prob_stay <- predict(mis_fit, newdata = pred_data, type = "response")
pred_data <- pred_data %>%
  mutate(drug = factor(drug, levels=c(0,1), labels=c("ddC","ddI")),
         prevOI = factor(prevOI, levels=c(0,1), labels=c("No Previous OI","Previous OI")))


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

# Save to Desktop
ggsave(
  filename = "C:/Users/Lenovo/Desktop/plot_5_Dropout_Probability.pdf",
  plot = p5,
  width = 9, height = 6, units = "in"
)










# --- 6. Model Diagnostic Plots ---
# Create a long-format data frame for diagnostic plots
res_long <- bind_rows(
  transmute(Data_complete, model = "GEE",  fitted = mu_gee,  resid = res_gee),
  transmute(Data_complete, model = "WGEE", fitted = mu_wgee, resid = res_wgee)
) %>%
  filter(is.finite(fitted), is.finite(resid)) %>%
  mutate(model = factor(model, levels = c("GEE","WGEE")))

# (Plot 6) Residuals vs. Fitted Plot
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

