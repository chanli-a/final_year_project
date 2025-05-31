library(ggplot2)
library(dplyr)

# 1. Read & prep data
df <- read.csv("output/metrics_final_coarse.csv")

df <- df[is.finite(df$ErrEqAb), ]
bins_K   <- sort(unique(df$T_bin_K))
labels_C <- c("10–12", "14–16", "18–20", "22–24",
              "26–28", "30–32", "34–36", "38–40")
df$T_bin_K <- factor(df$T_bin_K, levels = bins_K, labels = labels_C)

# 2. Define box width (same as used in geom_boxplot)
box_width <- 0.6

# 3. Build base plot without T-caps (not printed yet)
base_p <- ggplot(df, aes(x = T_bin_K, y = ErrEqAb)) +
  geom_boxplot(
    fill         = "lightgrey",
    colour       = "black",
    outlier.size = 1,
    width        = box_width,
    linewidth = 0.4
  ) +
  labs(x = "Temperature (°C)",
       y = "Equilibrium abundance deviation") +
  theme_classic() +
  theme(
    axis.text  = element_text(size = 12),
    axis.title = element_text(size = 14)
  )

# 4. Extract whisker endpoints and x positions from ggplot object
pb    <- ggplot_build(base_p)
stats <- pb$data[[1]]

# 5. Create data frame for T-caps (manually set width)
cap_df <- stats %>%
  mutate(
    x0 = x - box_width / 4,
    x1 = x + box_width / 4,
  )

# 6. Draw full-range plot with T-caps
p_full <- base_p +
  geom_segment(data = cap_df,
               aes(x = x0, xend = x1, y = ymin, yend = ymin),
               inherit.aes = FALSE,
               linewidth = 0.5) +
  geom_segment(data = cap_df,
               aes(x = x0, xend = x1, y = ymax, yend = ymax),
               inherit.aes = FALSE,
               linewidth = 0.5) +
  ggtitle("ErrEqAb by 4 °C Bin (full range)")

print(p_full)

# 7. Zoomed-in version
ylim_range <- range(c(stats$ymin, stats$ymax), na.rm = TRUE)
p_zoom <- p_full +
  coord_cartesian(ylim = ylim_range) +
  ggtitle("ErrEqAb by 4 °C Bin (zoomed to whiskers)")

print(p_zoom)









###############################################################################

####### reusable function for any metric (makes regular + zoomed in) ##########

###############################################################################



plot_temp_metric_boxplot <- function(df, temp_col, metric_col, 
                                     temp_labels = c("10–12", "14–16", "18–20", "22–24",
                                                     "26–28", "30–32", "34–36", "38–40"),
                                     box_width = 0.6) {
  
  # Ensure temperature column is a factor with labels
  bins <- sort(unique(df[[temp_col]]))
  df[[temp_col]] <- factor(df[[temp_col]], levels = bins, labels = temp_labels)
  
  # Filter for finite values in the metric
  df <- df[is.finite(df[[metric_col]]), ]
  
  # Base plot
  base_p <- ggplot(df, aes_string(x = temp_col, y = metric_col)) +
    geom_boxplot(
      fill         = "lightgrey",
      colour       = "black",
      outlier.size = 1,
      width        = box_width,
      linewidth = 0.4
    ) +
    labs(x = "Temperature (°C)",
         y = metric_col) +
    theme_classic() +
    theme(
      axis.text  = element_text(size = 12),
      axis.title = element_text(size = 14)
    )
  
  # Extract whiskers info
  pb    <- ggplot_build(base_p)
  stats <- pb$data[[1]] %>%
    mutate(
      x0 = x - box_width / 4,
      x1 = x + box_width / 4,
    )
  
  # Full range plot with T-caps
  p_full <- base_p +
    geom_segment(data = stats,
                 aes(x = x0, xend = x1, y = ymin, yend = ymin),
                 inherit.aes = FALSE,
                 linewidth = 0.5) +
    geom_segment(data = stats,
                 aes(x = x0, xend = x1, y = ymax, yend = ymax),
                 inherit.aes = FALSE,
                 linewidth = 0.5) +
    ggtitle(paste(metric_col, "by 4 °C Bin (full range)"))
  
  print(p_full)
  
  # Zoomed-in version to whiskers
  ylim_range <- range(c(stats$ymin, stats$ymax), na.rm = TRUE)
  p_zoom <- p_full +
    coord_cartesian(ylim = ylim_range) +
    ggtitle(paste(metric_col, "by 4 °C Bin (zoomed to whiskers)"))
  
  print(p_zoom)
}


plot_temp_metric_boxplot(df, temp_col = "T_bin_K", metric_col = "ErrEqAb")










###############################################################################

####### reusable function for any metric (makes regular / zoomed in - CHOOSE) ##########

###############################################################################









plot_temp_metric_boxplot_choose <- function(df, temp_col, metric_col, 
                                     temp_labels = c("10–12", "14–16", "18–20", "22–24",
                                                     "26–28", "30–32", "34–36", "38–40"),
                                     box_width = 0.6,
                                     zoom_only = FALSE,
                                     y_label = NULL) {
  
  # Use metric_col as default y-axis label if y_label not provided
  if (is.null(y_label)) {
    y_label <- metric_col
  }
  
  # Ensure temperature column is a factor with proper labels
  bins <- sort(unique(df[[temp_col]]))
  df[[temp_col]] <- factor(df[[temp_col]], levels = bins, labels = temp_labels)
  
  # Filter for finite values of the metric
  df <- df[is.finite(df[[metric_col]]), ]
  
  # Base plot
  base_p <- ggplot(df, aes_string(x = temp_col, y = metric_col)) +
    geom_boxplot(
      fill         = "lightgrey",
      colour       = "black",
      outlier.size = 1,
      width        = box_width,
      linewidth    = 0.4
    ) +
    labs(x = "Temperature (°C)", y = y_label) +
    theme_classic() +
    theme(
      axis.text  = element_text(size = 15),
      axis.title = element_text(size = 24)
    )
  
  # Extract whiskers info
  pb    <- ggplot_build(base_p)
  stats <- pb$data[[1]] %>%
    mutate(
      x0 = x - box_width / 4,
      x1 = x + box_width / 4,
    )
  
  # Base plot with T-caps
  p <- base_p +
    geom_segment(data = stats,
                 aes(x = x0, xend = x1, y = ymin, yend = ymin),
                 inherit.aes = FALSE,
                 linewidth = 0.5) +
    geom_segment(data = stats,
                 aes(x = x0, xend = x1, y = ymax, yend = ymax),
                 inherit.aes = FALSE,
                 linewidth = 0.5)
  
  # Add appropriate title
  title_suffix <- if (zoom_only) {
    " (zoomed to whiskers)"
  } else {
    " (full range)"
  }
  
  p <- p + ggtitle(paste(metric_col, "by 4 °C Bin", title_suffix))
  
  # Apply zoom if requested
  if (zoom_only) {
    ylim_range <- range(c(stats$ymin, stats$ymax), na.rm = TRUE)
    p <- p + coord_cartesian(ylim = ylim_range)
  }
  
  return(p)
}


# diversity 

plot_temp_metric_boxplot_choose(df, "T_bin_K", "jaccard",y_label = "Jaccard index")
plot_temp_metric_boxplot_choose(df, "T_bin_K", "bray_curtis",y_label = "Bray-Curtis dissimilarity")
plot_temp_metric_boxplot_choose(df, "T_bin_K", "bray_curtis",y_label = "Bray-Curtis index",zoom_only=TRUE)
plot_temp_metric_boxplot_choose(df, "T_bin_K", "shannon_mi",y_label = "MiCRM Shannon diversity index")
plot_temp_metric_boxplot_choose(df, "T_bin_K", "shannon_lv",y_label = "ELVM Shannon diversity index")


# equilibrium / trajectory abundance 

plot_temp_metric_boxplot_choose(df, "T_bin_K", "ErrEqAb",y_label = "Equilibrium abundance deviation log-ratio")
plot_temp_metric_boxplot_choose(df, "T_bin_K", "ErrEqAb",y_label = "Equilibrium abundance deviation log-ratio",zoom_only = TRUE)
plot_temp_metric_boxplot_choose(df, "T_bin_K", "ErrTraj",y_label = "Trajectory deviation log-ratio")
plot_temp_metric_boxplot_choose(df, "T_bin_K", "ErrTraj",y_label = "Trajectory deviation log-ratio", zoom_only=TRUE)


# stability / reactivity 

plot_temp_metric_boxplot_choose(df, "T_bin_K", "abs_stab_err",y_label = expression("|"*Re(lambda[dom](J[LV]))-Re(lambda[dom](J[Mi]))*"|"))
plot_temp_metric_boxplot_choose(df, "T_bin_K", "abs_react_err",y_label = expression("|"*Re(lambda[dom](H[LV]))-Re(lambda[dom](H[Mi]))*"|"))
# other plots are not box plots (continuous GLV vs MiCRM stability/reactivity) 


# timescale separation 
plot_temp_metric_boxplot_choose(df, "T_bin_K", "epsilon",y_label = expression(epsilon))
plot_temp_metric_boxplot_choose(df, "T_bin_K", "log10_eps_t_eq",y_label = expression(log[10](epsilon / t[eq])))


# hessian 
plot_temp_metric_boxplot_choose(df, "T_bin_K", "hessian_norm", y_label = "Hessian norm")


# non-normality
plot_temp_metric_boxplot_choose(df, "T_bin_K", "non_normality", y_label = expression(paste("|| ", J * J^{T} - J^{T} * J, " ||")))


# test ELVM stability / reactivity just to see what it looks like
plot_temp_metric_boxplot_choose(df, "T_bin_K", "stab_glv",y_label = expression(paste("ELVM ", Re(lambda[dom](J)))))
plot_temp_metric_boxplot_choose(df, "T_bin_K", "stab_glv",y_label = expression(paste("ELVM ", Re(lambda[dom](J)))), zoom_only=TRUE)
plot_temp_metric_boxplot_choose(df, "T_bin_K", "react_glv", 
                                y_label =expression(Re(lambda[dom](H[LV])))) +
  geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.4)

# box plots for stab/react spread (individual) 

plot_temp_metric_boxplot_choose(df, "T_bin_K", "stab_mi",y_label = expression(Re(lambda[dom](J[Mi]))))
plot_temp_metric_boxplot_choose(df, "T_bin_K", "stab_glv",y_label = expression(Re(lambda[dom](J[LV]))))
plot_temp_metric_boxplot_choose(df, "T_bin_K", "react_mi",y_label = expression(Re(lambda[dom](H[Mi]))))
plot_temp_metric_boxplot_choose(df, "T_bin_K", "react_glv",y_label = expression(Re(lambda[dom](H[LV]))))



#####################################################################################
############# NOW FOR THE LINE PLOTS ################################################
#####################################################################################




df_fine <- read_csv("output/metrics_final.csv") 




# timescale separation tauR vs tauC 


df_tau <- df_fine %>%
  group_by(T_C) %>%
  summarise(
    medC  = median(tau_C),
    medR  = median(tau_R)
  )

ggplot(df_tau, aes(x = T_C)) +
  # Median lines and points only
  geom_line(aes(y = medC, color = "Consumer median"), size = 0.8) +
  geom_line(aes(y = medR, color = "Resource median"), size = 0.8) +
  geom_point(aes(y = medC, color = "Consumer median"), size = 2) +
  geom_point(aes(y = medR, color = "Resource median"), size = 2) +
  scale_color_manual(
    name = NULL,
    values = c(
      "Consumer median" = "#005288",
      "Resource median" = "darkorange"   
    )
  ) +
  theme_classic() +
  labs(
    title = "Characteristic timescales of consumers and resources",
    x     = "Temperature (°C)",
    y     = "Characteristic timescales"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.title = element_text(size = 20),
    axis.text  = element_text(size = 15),
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 12),
    legend.position = c(0.75, 0.8),
    legend.background = element_rect(fill = alpha("white", 0.7), color = NA)
  )



# timescale separation epsilon 



df_eps <- df %>%
  group_by(T_C) %>%
  summarise(
    med = median(epsilon)
  )

ggplot(df_eps, aes(x = T_C, y = med)) +
  geom_line(color = "black", size = 0.8) +
  geom_point(color = "black", size = 2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  theme_classic() +
  labs(
    x = "Temperature (°C)",
    y = expression(epsilon)
  ) +
  theme(
    axis.title = element_text(size = 20),
    axis.text  = element_text(size = 15),
    legend.position = "none",
    plot.title = element_blank()
  )



# stab/react: plot ELVM and MiCRM ranges separately (probably supplementary) 


plot_overlay_single <- function(df, col, lab, ylabel, hlines = NULL, hline = NULL) {
  
  # Labels for legend
  line_lab <- paste0(lab, " median")
  fill_lab <- paste0(lab, " 90% range")
  
  # Summarise stats
  df_stats <- df %>%
    group_by(T_C) %>%
    summarise(
      med   = median(.data[[col]], na.rm = TRUE),
      low90 = quantile(.data[[col]], 0.05, na.rm = TRUE),
      high90 = quantile(.data[[col]], 0.95, na.rm = TRUE)
    )
  
  # Plot
  ggplot(df_stats, aes(x = T_C)) +
    # Ribbon with legend
    geom_ribbon(aes(ymin = low90, ymax = high90, fill = fill_lab), alpha = 0.4) +
    
    # Median line and points with legend
    geom_line(aes(y = med, color = line_lab), size = 0.8) +
    geom_point(aes(y = med, color = line_lab), size = 2, shape = ifelse(lab == "MiCRM", 16, 17)) +
    
    # Optional horizontal lines
    { if (!is.null(hline)) geom_hline(yintercept = hline, linetype = "dashed", color = "grey40") } +
    { if (!is.null(hlines)) lapply(hlines, function(h) geom_hline(yintercept = h, linetype = "dashed", color = "grey40")) } +
    
    # Manual legend colors (black for line, grey for ribbon)
    scale_color_manual(values = setNames("black", line_lab)) +
    scale_fill_manual(values = setNames("grey80", fill_lab)) +
    
    labs(
      x = "Temperature (°C)",
      y = ylabel,
      color = NULL,
      fill = NULL
    ) +
    
    theme_classic() +
    theme(
      axis.title = element_text(size = 24),
      axis.text  = element_text(size = 15),
      legend.text = element_text(size = 15),
      legend.position = c(0.2, 0.2),
      legend.direction = "vertical",
      legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
      legend.box = "vertical",
      legend.spacing.y = unit(5, "pt")
    )
}



# before plotting either MiCRM / ELVM, make sure to adjust legend position 
# for MiCRM plot, move legend to bottom left corner (0.2, 0.2) 
# for ELVM plot, move legend to top left corner (0.15, 0.85)


# single plots for stability (to show individual spread) 
plot_overlay_single(df, "stab_mi", "MiCRM", expression(Re(lambda[dom](J[Mi]))))
plot_overlay_single(df, "stab_glv", "ELVM", expression(Re(lambda[dom](J[LV]))))


# single plots for reactivity (to show individual spread) 
plot_overlay_single(df, "react_mi", "MiCRM", expression(Re(lambda[dom](H[Mi]))))
plot_overlay_single(df, "react_glv", "ELVM", expression(Re(lambda[dom](H[LV]))))


# stab/react, no shading, overlay 

plot_overlay <- function(df, col1, col2, lab1, lab2, ylabel, hlines = NULL, hline = NULL) {
  
  # Legend labels (median only now)
  line_lab1 <- paste0(lab1, " median")
  line_lab2 <- paste0(lab2, " median")
  
  # Colors (still distinguishable)
  color_vals <- setNames(c("#005288", "darkorange"), c(line_lab1, line_lab2))
  
  # Summarised data
  df_stats <- df %>%
    group_by(T_C) %>%
    summarise(
      med1 = median(.data[[col1]], na.rm = TRUE),
      med2 = median(.data[[col2]], na.rm = TRUE)
    )
  
  # Plot
  ggplot(df_stats, aes(x = T_C)) +
    # Median lines and points
    geom_line(aes(y = med1, color = line_lab1), size = 0.8) +
    geom_point(aes(y = med1, color = line_lab1), shape = 16, size = 2) +
    geom_line(aes(y = med2, color = line_lab2), size = 0.8) +
    geom_point(aes(y = med2, color = line_lab2), shape = 17, size = 2) +
    
    # Optional horizontal lines
    { if (!is.null(hline)) geom_hline(yintercept = hline, linetype = "dashed", color = "grey40") } +
    { if (!is.null(hlines)) lapply(hlines, function(h) geom_hline(yintercept = h, linetype = "dashed", color = "grey40")) } +
    
    # Manual legends
    scale_color_manual(values = color_vals) +
    
    labs(
      x = "Temperature (°C)",
      y = ylabel,
      color = NULL
    ) +
    
    theme_classic() +
    theme(
      axis.title = element_text(size = 20),
      axis.text  = element_text(size = 15),
      legend.text = element_text(size = 20),
      legend.position = c(0.2, 0.2),
      legend.direction = "vertical",
      legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
      legend.box = "vertical",
      legend.spacing.y = unit(5, "pt")
    )
}


plot_overlay(
  df = df,
  col1 = "stab_mi",
  col2 = "stab_glv",
  lab1 = "MiCRM",
  lab2 = "ELVM",
  ylabel = expression(Re(lambda[dom](J)))
)


plot_overlay(
  df = df,
  col1 = "react_mi",
  col2 = "react_glv",
  lab1 = "MiCRM",
  lab2 = "ELVM",
  ylabel = expression(Re(lambda[dom](H)))
)



# other metrics (not temp on x-axis): hessian/non-normality vs stability/reactivity 



# A reusable function to draw one scatter
plot_scatter <- function(df, xcol, ycol, xl, yl, title = "", log_x = FALSE, log_y = FALSE) {
  p <- ggplot(df, aes_string(x = xcol, y = ycol, color = "T_C")) +
    geom_point(alpha = 0.6, size = 2, shape = 16) +  # semi-transparent filled dots
    scale_color_gradientn(
      colours = c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"),
      name = "Temperature (°C)",
      limits = c(10, 40),
    ) +
    labs(
      x = xl,
      y = yl,
      title = title
    ) +
    theme_classic() +
    theme(
      plot.title   = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title   = element_text(size = 20),
      axis.text    = element_text(size = 15),
      legend.title = element_text(size = 15),
      legend.text  = element_text(size = 15)
    )
  
  # Apply log scales if requested
  if (log_x) p <- p + scale_x_log10()
  if (log_y) p <- p + scale_y_log10()
  
  return(p)
}




# Plot 1 - hermitian vs stability 
p1 <- plot_scatter(
  df,
  xcol  = "hessian_norm",
  ycol  = "abs_stab_err",
  xl    = "Log(Hessian norm)",
  yl    = expression("Log( |"*Re(lambda[dom](J[LV]))-Re(lambda[dom](J[Mi]))*"| )"),
  log_x=TRUE,
  log_y=TRUE
)

p1 <- p1 + geom_smooth(method = "lm", se = FALSE, color = "black", size=0.7)

model1 <- lm(log10(abs_stab_err) ~ log10(hessian_norm), data = df)
summ1 <- summary(model1)
r2_1 <- summ1$r.squared
p_slope1 <- summ1$coefficients["log10(hessian_norm)", "Pr(>|t|)"]

p1 <- p1 + annotate(
  geom  = "text",
  x     = Inf,                # place at right edge
  y     = Inf,                # place at top edge
  label = bquote(R^2 == .(format(r2_1, digits = 2))),
  hjust = 1.1,                # nudge inward horizontally
  vjust = 1.5,                # nudge inward vertically
  size  = 5                   # adjust text size
)

print(p1)
cat(sprintf("R² = %.3f\np-value (slope) = %.3g\n", r2_1, p_slope1))




# Plot 2 - hermitian vs reactivity 
p2 <- plot_scatter(
  df,
  xcol  = "hessian_norm",
  ycol  = "abs_react_err",
  xl    = "Log(Hessian norm)",
  yl    = expression("Log( |"*Re(lambda[dom](H[LV]))-Re(lambda[dom](H[Mi]))*"| )"),
  log_x=TRUE,
  log_y=TRUE
)
print(p2)

# Plot 3 - non normality vs stability 
p3 <- plot_scatter(
  df,
  xcol  = "non_normality",
  ycol  = "abs_stab_err",
  xl    = expression(paste("Log(|| ", J * J^{T} - J^{T} * J, " ||)")),
  yl    = expression("Log( |"*Re(lambda[dom](J[LV]))-Re(lambda[dom](J[Mi]))*"| )"),
  log_x=TRUE,
  log_y=TRUE
)
print(p3)

# Plot 4 - non normality vs reactivity 
p4 <- plot_scatter(
  df,
  xcol  = "non_normality",
  ycol  = "abs_react_err",
  xl    = expression(paste("Log(|| ", J * J^{T} - J^{T} * J, " ||)")),
  yl    = expression("Log( |"*Re(lambda[dom](H[LV]))-Re(lambda[dom](H[Mi]))*"| )"),
  log_x=TRUE,
  log_y=TRUE
)
print(p4)




###########

# just out of curiosity 




# hessian vs epsilon (does HOI drive timescale sep?)
p5 <- plot_scatter(
  df,
  xcol  = "hessian_norm",
  ycol  = "epsilon",
  xl    = "Log(Hessian norm)",
  yl    = expression(epsilon),
  log_x = TRUE
)

# Add horizontal line and set y-axis breaks at 1-unit intervals
p5 <- p5 +
  geom_hline(yintercept = 1, linetype = "dotted", color = "grey40") +
  scale_y_continuous(breaks = seq(floor(min(df$epsilon, na.rm = TRUE)),
                                  ceiling(max(df$epsilon, na.rm = TRUE)), 
                                  by = 1))

print(p5)


# non normality vs epsilon (does interconnectedness drive timescale sep?)
p6 <- plot_scatter(
  df,
  xcol  = "non_normality",
  ycol  = "epsilon",
  xl    = expression(paste("Log(|| ", J * J^{T} - J^{T} * J, " ||)")),
  yl    = expression(epsilon),
  log_x=TRUE
)

p6 <- p6 +
  geom_hline(yintercept = 1, linetype = "dotted", color = "grey40") +
  scale_y_continuous(breaks = seq(floor(min(df$epsilon, na.rm = TRUE)),
                                  ceiling(max(df$epsilon, na.rm = TRUE)), 
                                  by = 1))

print(p6)







