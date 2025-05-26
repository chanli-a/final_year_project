install.packages("tidyverse")

library(tidyverse)
library(ggplot2)

setwd("C:/Users/Chan Li/OneDrive - Imperial College London/Third Year/FYP/Code-related/final_year_project/temp_abundance")
getwd()

# ─── 1) Load & preprocess ───────────────────────────────────────────────────
df <- read_csv("output/metrics_N50.csv") %>%
  mutate(T_C = T_K - 273.15)

df_traj <- read_csv("output/traj_N50.csv") %>%
  mutate(T_C = T_K - 273.15)


# Equilibrium abundance error

ggplot(df_err, aes(x = T_C)) +
  geom_ribbon(aes(ymin = low95, ymax = high95, fill = "95% interval"), alpha = 0.3) +
  geom_ribbon(aes(ymin = p25,  ymax = p75,  fill = "50% IQR"),       alpha = 0.6) +
  geom_line(aes(y = med, color = "Median"),      size = 0.8) +
  geom_point(aes(y = med, color = "Median", shape = "Median"), size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(
    name   = "",
    values = c("95% interval" = "lightblue", "50% IQR" = "#c49ec4"),
    guide  = guide_legend(order = 2)
  ) +
  scale_color_manual(
    name   = "",
    values = c("Median" = "black"),
    guide  = guide_legend(order = 1)
  ) +
  scale_shape_manual(
    name   = "",
    values = c("Median" = 16),
    guide  = guide_legend(order = 1)
  ) +
  theme_classic() +
  labs(x = "Temperature (°C)", y = "Error") +
  theme(
    axis.title      = element_text(size = 16),
    axis.text       = element_text(size = 14),
    legend.title    = element_text(size = 14),
    legend.text     = element_text(size = 12),
    legend.position = c(0.2, 0.2),
    legend.background = element_rect(fill = alpha("white", 0.7), color = NA)
  )

# updated equilibrium abundance, without 95% interval 

ggplot(df_err, aes(x = T_C)) +
  geom_ribbon(aes(ymin = p25, ymax = p75, fill = "50% IQR"), alpha = 0.6) +
  geom_line(aes(y = med, color = "Median", group = 1), size = 0.8) +
  geom_point(aes(y = med, color = "Median", shape = "Median"), size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(
    name = "",
    values = c("50% IQR" = "#c49ec4"),
    guide = guide_legend(order = 2)
  ) +
  scale_color_manual(
    name = "",
    values = c("Median" = "black"),
    guide = guide_legend(order = 1)
  ) +
  scale_shape_manual(
    name = "",
    values = c("Median" = 16),
    guide = guide_legend(order = 1)
  ) +
  theme_classic() +
  labs(x = "Temperature (°C)", y = "Error") +
  theme(
    axis.title        = element_text(size = 16),
    axis.text         = element_text(size = 14),
    legend.title      = element_text(size = 14),
    legend.text       = element_text(size = 12),
    legend.position   = c(0.2, 0.2),
    legend.spacing.y  = unit(0.0, 'cm'),
    legend.background = element_rect(fill = alpha("white", 0.7), color = NA)
  ) +
  guides(
    fill = guide_legend(order = 2),
    color = guide_legend(order = 1),
    shape = guide_legend(order = 1)
  )



# trajectory deviation

df_td <- df_traj %>%
  group_by(T_C) %>%
  summarise(
    med    = median(ErrTraj),
    p25    = quantile(ErrTraj, 0.25),
    p75    = quantile(ErrTraj, 0.75),
    low95  = quantile(ErrTraj, 0.025),
    high95 = quantile(ErrTraj, 0.975)
  )

ggplot(df_td, aes(x = T_C)) +
  geom_ribbon(aes(ymin = low95, ymax = high95, fill = "95% interval"), alpha = 0.3) +
  geom_ribbon(aes(ymin = p25, ymax = p75, fill = "50% IQR"), alpha = 0.6) +
  geom_line(aes(y = med, color = "Median"), size = 0.8) +
  geom_point(aes(y = med, color = "Median", shape = "Median"), size = 2) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  scale_fill_manual(name = "", values = c("95% interval" = "lightblue", "50% IQR" = "#c49ec4")) +
  scale_color_manual(name = "", values = c("Median" = "black")) +
  scale_shape_manual(name = "", values = c("Median" = 16)) +
  theme_classic() +
  labs(x = "Temperature (°C)", y = "ErrTraj") +
  theme(
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 12),
    legend.position = c(0.2, 0.2),
    legend.background = element_rect(fill = alpha("white", 0.7), color = NA)
  )




# species overlaps 
df_ov <- df %>%
  group_by(T_C) %>%
  summarise(
    mean   = mean(overlap,      na.rm = TRUE),
    p25    = quantile(overlap, 0.25, na.rm = TRUE),
    p75    = quantile(overlap, 0.75, na.rm = TRUE),
    low95  = quantile(overlap, 0.025, na.rm = TRUE),
    high95 = quantile(overlap, 0.975, na.rm = TRUE)
  ) %>% ungroup()

# Outliers
df_ov_out <- df %>%
  group_by(T_C) %>%
  filter(overlap < quantile(overlap, 0.025, na.rm=TRUE) |
           overlap > quantile(overlap, 0.975, na.rm=TRUE)) %>%
  ungroup()

# Plot
p2 <- ggplot() +
  geom_ribbon(data = df_ov,
              aes(x = T_C, ymin = low95, ymax = high95, fill = "95% interval"),
              alpha = 0.3) +
  geom_ribbon(data = df_ov,
              aes(x = T_C, ymin = p25, ymax = p75, fill = "50% IQR"),
              alpha = 0.6) +
  geom_line(data = df_ov,
            aes(x = T_C, y = mean, color = "Mean"),
            size = 0.8) +
  geom_point(data = df_ov,
             aes(x = T_C, y = mean, shape = "Mean", color = "Mean"),
             size = 2) +
  geom_point(data = df_ov_out,
             aes(x = T_C, y = overlap, shape = "Outliers"),
             color = "lightblue", size = 2, alpha = 0.7) +
  scale_fill_manual(name = "", values = c("50% IQR" = "#c49ec4", "95% interval" = "lightblue")) +
  scale_color_manual(name = "", values = c("Mean" = "black")) +
  scale_shape_manual(name = "", values = c("Mean" = 16, "Outliers" = 17)) +
  theme_classic() +
  labs(
    title = "Species Overlap",
    x     = "Temperature (°C)",
    y     = "Number of overlapping species"
  )
print(p2)



# jaccard index 

# Update summary dataframe to include median
df_j <- df %>%
  group_by(T_C) %>%
  summarise(
    median = median(jaccard,     na.rm = TRUE),
    p25    = quantile(jaccard, 0.25, na.rm = TRUE),
    p75    = quantile(jaccard, 0.75, na.rm = TRUE),
    low95  = quantile(jaccard, 0.025, na.rm = TRUE),
    high95 = quantile(jaccard, 0.975, na.rm = TRUE)
  ) %>% ungroup()

# Plot without outliers, using median
p3 <- ggplot() +
  geom_ribbon(data = df_j,
              aes(x = T_C, ymin = low95, ymax = high95, fill = "95% interval"),
              alpha = 0.3) +
  geom_ribbon(data = df_j,
              aes(x = T_C, ymin = p25, ymax = p75, fill = "50% IQR"),
              alpha = 0.6) +
  geom_line(data = df_j,
            aes(x = T_C, y = median, color = "Median"),
            size = 0.8) +
  geom_point(data = df_j,
             aes(x = T_C, y = median, shape = "Median", color = "Median"),
             size = 2) +
  scale_fill_manual(
    name = "",
    values = c("50% IQR" = "#c49ec4", "95% interval" = "lightblue")
  ) +
  scale_color_manual(
    name = "",
    values = c("Median" = "black")
  ) +
  scale_shape_manual(
    name = "",
    values = c("Median" = 16)
  ) +
  geom_hline(yintercept = 1.0, linetype = "dashed") +
  coord_cartesian(ylim = c(NA, 1.2)) +  # Limit y-axis max to 1.2
  theme_classic() +
  labs(
    title = "Jaccard Similarity",
    x     = "Temperature (°C)",
    y     = "Jaccard Index"
  )

print(p3)




df_bc <- df %>%
  group_by(T_C) %>%
  summarise(
    mean   = mean(bray_curtis),
    p25    = quantile(bray_curtis, 0.25),
    p75    = quantile(bray_curtis, 0.75),
    low95  = quantile(bray_curtis, 0.025),
    high95 = quantile(bray_curtis, 0.975)
  )

ggplot(df_bc, aes(x = T_C)) +
  geom_ribbon(aes(ymin = low95, ymax = high95),
              fill = "lightblue", alpha = 0.3) +
  geom_ribbon(aes(ymin = p25,  ymax = p75),
              fill = "#c49ec4", alpha = 0.6) +
  geom_line(aes(y = mean), size = 0.8) +
  geom_point(aes(y = mean), size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_classic() +
  labs(
    title = "Bray–Curtis Dissimilarity",
    x     = "Temperature (°C)",
    y     = "Bray–Curtis"
  )


df_sh <- df %>%
  group_by(T_C) %>%
  summarise(
    m_mi    = mean(shannon_mi),
    p25_mi  = quantile(shannon_mi, 0.25),
    p75_mi  = quantile(shannon_mi, 0.75),
    low95_mi= quantile(shannon_mi, 0.025),
    high95_mi=quantile(shannon_mi, 0.975),
    m_lv    = mean(shannon_lv),
    p25_lv  = quantile(shannon_lv, 0.25),
    p75_lv  = quantile(shannon_lv, 0.75),
    low95_lv= quantile(shannon_lv, 0.025),
    high95_lv=quantile(shannon_lv, 0.975)
  )

ggplot(df_sh, aes(x = T_C)) +
  # MiCRM band
  geom_ribbon(aes(ymin = low95_mi, ymax = high95_mi),
              fill = "lightblue", alpha = 0.3) +
  geom_ribbon(aes(ymin = p25_mi,  ymax = p75_mi),
              fill = "#00bfff", alpha = 0.6) +
  # GLV band
  geom_ribbon(aes(ymin = low95_lv, ymax = high95_lv),
              fill = "lightblue", alpha = 0.3) +
  geom_ribbon(aes(ymin = p25_lv,  ymax = p75_lv),
              fill = "#ff8c00", alpha = 0.6) +
  # Lines
  geom_line(aes(y = m_mi), color = "#00bfff", size = 0.8) +
  geom_line(aes(y = m_lv), color = "#ff8c00", size = 0.8) +
  geom_point(aes(y = m_mi), color = "#00bfff", size = 2) +
  geom_point(aes(y = m_lv), color = "#ff8c00", size = 2) +
  theme_classic() +
  labs(
    title = "Shannon Diversity (MiCRM vs. GLV)",
    x     = "Temperature (°C)",
    y     = "Shannon Index",
    color = ""
  )


# Stability overlay
df_stab <- df %>%
  group_by(T_C) %>%
  summarise(
    m_mi    = mean(stab_mi),
    p25_mi  = quantile(stab_mi, 0.25),
    p75_mi  = quantile(stab_mi, 0.75),
    low95_mi= quantile(stab_mi, 0.025),
    high95_mi=quantile(stab_mi, 0.975),
    m_lv    = mean(stab_glv),
    p25_lv  = quantile(stab_glv, 0.25),
    p75_lv  = quantile(stab_glv, 0.75),
    low95_lv= quantile(stab_glv, 0.025),
    high95_lv=quantile(stab_glv, 0.975)
  )

ggplot(df_stab, aes(x = T_C)) +
  geom_ribbon(aes(ymin = low95_mi, ymax = high95_mi),
              fill = "#87cefa", alpha = 0.3) +
  geom_ribbon(aes(ymin = p25_mi,  ymax = p75_mi),
              fill = "#00bfff", alpha = 0.6) +
  geom_ribbon(aes(ymin = low95_lv, ymax = high95_lv),
              fill = "#ffe4b5", alpha = 0.3) +
  geom_ribbon(aes(ymin = p25_lv,  ymax = p75_lv),
              fill = "#ff8c00", alpha = 0.6) +
  geom_line(aes(y = m_mi), color = "#00bfff", size = 0.8) +
  geom_line(aes(y = m_lv), color = "#ff8c00", size = 0.8) +
  geom_point(aes(y = m_mi), color = "#00bfff", size = 2) +
  geom_point(aes(y = m_lv), color = "#ff8c00", size = 2) +
  theme_classic() +
  labs(
    title = "Stability (MiCRM vs. GLV)",
    x     = "Temperature (°C)",
    y     = "Leading Eigenvalue"
  )

# Absolute Stability Error
df_absstab <- df %>%
  group_by(T_C) %>%
  summarise(
    mean   = mean(abs_stab_err),
    p25    = quantile(abs_stab_err, 0.25),
    p75    = quantile(abs_stab_err, 0.75),
    low95  = quantile(abs_stab_err, 0.025),
    high95 = quantile(abs_stab_err, 0.975)
  )

ggplot(df_absstab, aes(x = T_C)) +
  geom_ribbon(aes(ymin = low95, ymax = high95),
              fill = "lightblue", alpha = 0.3) +
  geom_ribbon(aes(ymin = p25,  ymax = p75),
              fill = "#c49ec4", alpha = 0.6) +
  geom_line(aes(y = mean), size = 0.8) +
  geom_point(aes(y = mean), size = 2) +
  theme_classic() +
  labs(
    title = "Absolute Stability Error",
    x     = "Temperature (°C)",
    y     = "Error"
  )


# Reactivity overlay
df_react <- df %>%
  group_by(T_C) %>%
  summarise(
    m_mi    = mean(react_mi),
    p25_mi  = quantile(react_mi, 0.25),
    p75_mi  = quantile(react_mi, 0.75),
    low95_mi= quantile(react_mi, 0.025),
    high95_mi=quantile(react_mi, 0.975),
    m_lv    = mean(react_glv),
    p25_lv  = quantile(react_glv, 0.25),
    p75_lv  = quantile(react_glv, 0.75),
    low95_lv= quantile(react_glv, 0.025),
    high95_lv=quantile(react_glv, 0.975)
  )

ggplot(df_react, aes(x = T_C)) +
  geom_ribbon(aes(ymin = low95_mi, ymax = high95_mi),
              fill = "#87cefa", alpha = 0.3) +
  geom_ribbon(aes(ymin = p25_mi,  ymax = p75_mi),
              fill = "#00bfff", alpha = 0.6) +
  geom_ribbon(aes(ymin = low95_lv, ymax = high95_lv),
              fill = "#ffe4b5", alpha = 0.3) +
  geom_ribbon(aes(ymin = p25_lv,  ymax = p75_lv),
              fill = "#ff8c00", alpha = 0.6) +
  geom_line(aes(y = m_mi), color = "#00bfff", size = 0.8) +
  geom_line(aes(y = m_lv), color = "#ff8c00", size = 0.8) +
  geom_point(aes(y = m_mi), color = "#00bfff", size = 2) +
  geom_point(aes(y = m_lv), color = "#ff8c00", size = 2) +
  theme_classic() +
  labs(
    title = "Reactivity (MiCRM vs. GLV)",
    x     = "Temperature (°C)",
    y     = "Hermitian Leading Eigenvalue"
  )

# Absolute Reactivity Error
df_absreact <- df %>%
  group_by(T_C) %>%
  summarise(
    mean   = mean(abs_react_err),
    p25    = quantile(abs_react_err, 0.25),
    p75    = quantile(abs_react_err, 0.75),
    low95  = quantile(abs_react_err, 0.025),
    high95 = quantile(abs_react_err, 0.975)
  )

ggplot(df_absreact, aes(x = T_C)) +
  geom_ribbon(aes(ymin = low95, ymax = high95),
              fill = "lightblue", alpha = 0.3) +
  geom_ribbon(aes(ymin = p25,  ymax = p75),
              fill = "#c49ec4", alpha = 0.6) +
  geom_line(aes(y = mean), size = 0.8) +
  geom_point(aes(y = mean), size = 2) +
  theme_classic() +
  labs(
    title = "Absolute Reactivity Error",
    x     = "Temperature (°C)",
    y     = "Error"
  )


# Hessian norm
df_hess <- df %>%
  group_by(T_C) %>%
  summarise(
    mean   = mean(hessian_norm),
    p25    = quantile(hessian_norm, 0.25),
    p75    = quantile(hessian_norm, 0.75),
    low95  = quantile(hessian_norm, 0.025),
    high95 = quantile(hessian_norm, 0.975)
  )

ggplot(df_hess, aes(x = T_C)) +
  geom_ribbon(aes(ymin = low95, ymax = high95),
              fill = "lightblue", alpha = 0.3) +
  geom_ribbon(aes(ymin = p25,  ymax = p75),
              fill = "#c49ec4", alpha = 0.6) +
  geom_line(aes(y = mean), size = 0.8) +
  geom_point(aes(y = mean), size = 2) +
  theme_classic() +
  labs(
    title = "Hessian Frobenius Norm",
    x     = "Temperature (°C)",
    y     = "||H||"
  )

# Non-normality
df_nn <- df %>%
  group_by(T_C) %>%
  summarise(
    mean   = mean(non_normality),
    p25    = quantile(non_normality, 0.25),
    p75    = quantile(non_normality, 0.75),
    low95  = quantile(non_normality, 0.025),
    high95 = quantile(non_normality, 0.975)
  )

ggplot(df_nn, aes(x = T_C)) +
  geom_ribbon(aes(ymin = low95, ymax = high95),
              fill = "lightblue", alpha = 0.3) +
  geom_ribbon(aes(ymin = p25,  ymax = p75),
              fill = "#c49ec4", alpha = 0.6) +
  geom_line(aes(y = mean), size = 0.8) +
  geom_point(aes(y = mean), size = 2) +
  theme_classic() +
  labs(
    title = "Matrix Non-normality",
    x     = "Temperature (°C)",
    y     = "||[J, Jᵀ]||"
  )


# Epsilon
df_eps <- df %>%
  group_by(T_C) %>%
  summarise(
    med    = median(epsilon),
    p25    = quantile(epsilon, 0.25),
    p75    = quantile(epsilon, 0.75),
    low95  = quantile(epsilon, 0.025),
    high95 = quantile(epsilon, 0.975)
  )

ggplot(df_eps, aes(x = T_C)) +
  # 95% CI ribbon
  geom_ribbon(aes(ymin = low95, ymax = high95, fill = "95% CI"), alpha = 0.3) +
  # IQR ribbon
  geom_ribbon(aes(ymin = p25, ymax = p75, fill = "IQR"), alpha = 0.6) +
  # Median line and points
  geom_line(aes(y = med, color = "Median"), size = 0.8) +
  geom_point(aes(y = med, color = "Median"), size = 2) +
  # Reference line at 1
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  # Custom color/fill scales
  scale_fill_manual(
    name = "Shaded areas",
    values = c(
      "95% CI" = "lightblue",
      "IQR"    = "#c49ec4"
    )
  ) +
  scale_color_manual(
    name = "Line",
    values = c(
      "Median" = "black"  # dark blue for visibility
    )
  ) +
  # Theme and formatting
  theme_classic() +
  labs(
    title = "Consumer-resource timescale separation",
    x     = "Temperature (°C)",
    y     = expression(epsilon)
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 12),
    legend.position = c(0.75, 0.8),
    legend.background = element_rect(fill = alpha("white", 0.7), color = NA)
  )



# Log10(ε/t_eq)
df_loge <- df %>%
  group_by(T_C) %>%
  summarise(
    med    = median(log10_eps_t_eq),
    p25    = quantile(log10_eps_t_eq, 0.25),
    p75    = quantile(log10_eps_t_eq, 0.75),
    low95  = quantile(log10_eps_t_eq, 0.025),
    high95 = quantile(log10_eps_t_eq, 0.975)
  )

ggplot(df_loge, aes(x = T_C)) +
  # 95% CI ribbon
  geom_ribbon(aes(ymin = low95, ymax = high95, fill = "95% CI"), alpha = 0.3) +
  # IQR ribbon
  geom_ribbon(aes(ymin = p25, ymax = p75, fill = "IQR"), alpha = 0.6) +
  # Median line and points
  geom_line(aes(y = med, color = "Median"), size = 0.8) +
  geom_point(aes(y = med, color = "Median"), size = 2) +
  # Custom color/fill scales
  scale_fill_manual(
    name = "Shaded areas",
    values = c(
      "95% CI" = "lightblue",
      "IQR"    = "#c49ec4"
    )
  ) +
  scale_color_manual(
    name = "Line",
    values = c(
      "Median" = "black"  # dark blue
    )
  ) +
  # Theme and formatting
  theme_classic() +
  labs(
    title = "Consumer-resource timescale separation relative to system equilibration time",
    x     = "Temperature (°C)",
    y     = expression(log[10](epsilon / t[eq]))
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 12),
    legend.position = c(0.82, 0.8),
    legend.background = element_rect(fill = alpha("white", 0.7), color = NA)
  )




# timescale separation tauC and tauR 


df_tau <- df %>%
  group_by(T_C) %>%
  summarise(
    medC    = median(tau_C),
    p25C    = quantile(tau_C, 0.25),
    p75C    = quantile(tau_C, 0.75),
    low90C  = quantile(tau_C, 0.05),
    high90C = quantile(tau_C, 0.95),
    
    medR    = median(tau_R),
    p25R    = quantile(tau_R, 0.25),
    p75R    = quantile(tau_R, 0.75),
    low90R  = quantile(tau_R, 0.05),
    high90R = quantile(tau_R, 0.95)
  )

ggplot(df_tau, aes(x = T_C)) +
  # 90% ribbons
  geom_ribbon(aes(ymin = low90C, ymax = high90C, fill = "Consumer 90% CI"), alpha = 0.3) +
  geom_ribbon(aes(ymin = low90R, ymax = high90R, fill = "Resource 90% CI"), alpha = 0.3) +
  # 50% ribbons (IQR)
  geom_ribbon(aes(ymin = p25C, ymax = p75C, fill = "Consumer IQR"), alpha = 0.6) +
  geom_ribbon(aes(ymin = p25R, ymax = p75R, fill = "Resource IQR"), alpha = 0.6) +
  # Median lines
  geom_line(aes(y = medC, color = "Consumer Median"), size = 0.8) +
  geom_line(aes(y = medR, color = "Resource Median"), size = 0.8) +
  geom_point(aes(y = medC, color = "Consumer Median"), size = 2) +
  geom_point(aes(y = medR, color = "Resource Median"), size = 2) +
  scale_fill_manual(
    name = "Shaded Areas",
    values = c(
      "Consumer 90% CI" = "lightblue",
      "Consumer IQR"    = "#c49ec4",      # darker blue
      "Resource 90% CI" = "#ffdd99",      # lighter orange
      "Resource IQR"    = "#cc5500"       # warm orange/brown
    )
  ) +
  scale_color_manual(
    name = "Median Lines",
    values = c(
      "Consumer Median" = "#005288",   # dark steel blue
      "Resource Median" = "#8b4513"    # saddle brown
    )
  ) +
  theme_classic() +
  labs(
    title = "Characteristic Timescales of Consumers and Resources",
    x     = "Temperature (°C)",
    y     = "Characteristic Timescale"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 12)
  )

# timescale separation version with only 50% IQR

df_tau <- df %>%
  group_by(T_C) %>%
  summarise(
    medC  = median(tau_C),
    p25C  = quantile(tau_C, 0.25),
    p75C  = quantile(tau_C, 0.75),
    medR  = median(tau_R),
    p25R  = quantile(tau_R, 0.25),
    p75R  = quantile(tau_R, 0.75)
  )

ggplot(df_tau, aes(x = T_C)) +
  # 50% ribbons (IQR only)
  geom_ribbon(aes(ymin = p25C, ymax = p75C, fill = "Consumer IQR"), alpha = 0.6) +
  geom_ribbon(aes(ymin = p25R, ymax = p75R, fill = "Resource IQR"), alpha = 0.6) +
  # Median lines and points
  geom_line(aes(y = medC, color = "Consumer median"), size = 0.8) +
  geom_line(aes(y = medR, color = "Resource median"), size = 0.8) +
  geom_point(aes(y = medC, color = "Consumer median"), size = 2) +
  geom_point(aes(y = medR, color = "Resource median"), size = 2) +
  scale_fill_manual(
    name = "Shaded Areas",
    values = c(
      "Consumer IQR"  = "lightblue",   # muted purple
      "Resource IQR"  = "#cc5500"    # warm orange-brown
    )
  ) +
  scale_color_manual(
    name = "Lines",
    values = c(
      "Consumer median" = "#005288",  # dark steel blue
      "Resource median" = "#8b4513"   # saddle brown
    )
  ) +
  theme_classic() +
  labs(
    title = "Characteristic timescales of consumers and resources",
    x     = "Temperature (°C)",
    y     = "Characteristic timescale"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 12),
    legend.position = c(0.75, 0.8),
    legend.background = element_rect(fill = alpha("white", 0.7), color = NA)
  )



