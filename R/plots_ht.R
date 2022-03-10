my_theme <- theme_linedraw() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 14, face = "bold")
  )


## Cross-validation Plot
ht_cv_plot <- function(cv_df) {

  cv_plot <- cv_df %>%
    ggplot(aes(x = lambda, y = loglik)) +
    geom_line(size = 1) + geom_point(size = 2) +
    geom_vline(xintercept = best_lambda,
               size = 1.2, col = "red", linetype = "dashed") +
    geom_label(aes(x = best_lambda, y = max(loglik),
                   label = signif(best_lambda, 2)),
               colour = "white", fill = "red", fontface = "bold" ) +
    scale_y_continuous(breaks = number_ticks(8)) +
    scale_x_continuous(breaks = number_ticks(10)) +
    ggtitle(paste0(associates, "|", conditioning, " ", phase_name, " H&T Cross-Validation"),
            subtitle = "Roughness Penalty Selection") +
    labs(x = expression(lambda), y = "Negative Log-likelihood \n (log10 Scale)") +
    my_theme

  return(cv_plot)

}

ht_cv_plot_jointly <- function(cv_df, best_lambda, conditioning, phase_name, lead_id) {

  cv_plot <- cv_df %>%
    ggplot(aes(x = lambda, y = loglik)) +
    geom_line(size = 1) + geom_point(size = 2) +
    geom_vline(xintercept = best_lambda,
               size = 1.2, col = "red", linetype = "dashed") +
    geom_label(aes(x = best_lambda, y = max(loglik),
                   label = signif(best_lambda, 2)),
               colour = "white", fill = "red", fontface = "bold" ) +
    scale_y_continuous(breaks = number_ticks(8)) +
    scale_x_continuous(breaks = number_ticks(10)) +
    ggtitle(paste0("All Channels|", conditioning, " ", phase_name, " H&T Cross-Validation"),
            subtitle = paste0("Roughness Penalty Selection / Lag = ", lead_id)) +
    labs(x = expression(lambda), y = "Negative Log-likelihood \n (log10 Scale)") +
    my_theme

  return(cv_plot)

}


## Alpha Stability Plot
alpha_stb <- function(boot_summary, boot_nep, nep, tau, associates, phase_name, lead_id) {

  stb_brks <- c(seq(nep[1], tau, length.out = 4),
                seq(tau, nep[2], length.out = 4)[-1])

  blks_n <- length(unique(boot_summary$par)) - 3

  stb_df <- boot_summary %>%
    dplyr::filter(par %notin% c("beta", "mu", "sigma")) %>%
    mutate(par = factor(par, ordered = TRUE,
                        levels = paste0("alpha_", 1:blks_n),
                        labels = paste0("Block ", 1:blks_n))) %>%
    dplyr::select(all_of(associates), par) %>%
    cbind(boot_nep) %>%
    pivot_longer(cols = all_of(associates), names_to = "Channel", values_to = "Estimates")

  scaleFUN <- function(x) sprintf("%.3f", x)
  stb_df %>%
    dplyr::filter(Channel == associates) %>%
    ggplot(aes(x = boot_nep, y = Estimates)) +
    geom_point() +
    geom_smooth(method = "loess", formula = "y ~ x", col = "red") +
    scale_x_continuous(breaks = stb_brks, labels = scaleFUN, expand = c(0, 0)) +
    scale_y_continuous(breaks = number_ticks(5), expand = c(0, 0.05)) +
    facet_wrap(~par, nrow = 2) +
    ggtitle(bquote(bold("Channel"~.(associates)~"|T3"~.(phase_name)~"H&T"~alpha~"-stability by Threshold")),
            subtitle = paste0("500 Bootstrap Samples / Lag = ", lead_id)) +
    labs(y = expression(hat(alpha)), x = "NEP") +
    my_theme +
    theme(panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1)) -> alpha_plot

  return(alpha_plot)

}


## Beta Stability Plot
beta_stb <- function(boot_summary, boot_nep, nep, tau, associates, phase_name, lead_id) {

  stb_brks <- c(seq(nep[1], tau, length.out = 6),
                seq(tau, nep[2], length.out = 6)[-1])

  stb_df <- boot_summary %>%
    dplyr::filter(par == "beta") %>%
    dplyr::select(all_of(associates), par) %>%
    cbind(boot_nep) %>%
    pivot_longer(cols = all_of(associates), names_to = "Channel", values_to = "Estimates")

  scaleFUN <- function(x) sprintf("%.4f", x)
  stb_df %>%
    ggplot(aes(x = boot_nep, y = Estimates)) +
    geom_point() +
    geom_smooth(method = "loess", formula = "y ~ x", col = "red") +
    scale_x_continuous(breaks = stb_brks, labels = scaleFUN, expand = c(0, 0)) +
    scale_y_continuous(breaks = number_ticks(8), expand = c(0, 0.05)) +
    ggtitle(bquote(bold("Channel"~.(associates)~"|T3"~.(phase_name)~"H&T"~beta~"-stability by Threshold")),
            subtitle = paste0("500 Bootstrap Samples / Lag = ", lead_id)) +
    labs(y = expression(hat(beta)), x = "NEP") +
    my_theme +
    theme(panel.grid.minor = element_blank()) -> beta_plot

  return(beta_plot)

}


violin_plot <- function(boot_summary, boot_nep, nep, tau, associates, phase_name) {

  stb_brks <- c(seq(nep[1], tau, length.out = 6),
                seq(tau, nep[2], length.out = 6)[-1])

  stb_df <- boot_summary %>%
    dplyr::filter(par == "beta") %>%
    dplyr::select(all_of(associates), par) %>%
    cbind(boot_nep) %>%
    pivot_longer(cols = all_of(associates), names_to = "Channel", values_to = "Estimates")

  scaleFUN <- function(x) sprintf("%.4f", x)
  stb_df %>%
    ggplot(aes(x = boot_nep, y = Estimates)) +
    geom_point() +
    geom_smooth(method = "loess", formula = "y ~ x", col = "red") +
    scale_x_continuous(breaks = stb_brks, labels = scaleFUN, expand = c(0, 0)) +
    scale_y_continuous(breaks = number_ticks(8), expand = c(0, 0.05)) +
    ggtitle(bquote(bold("Channel"~.(associates)~.(phase_name)~"H&T"~beta~"-stability by Threshold")),
            subtitle = paste0("500 Bootstrap Samples")) +
    labs(y = expression(hat(beta)), x = "NEP") +
    my_theme +
    theme(panel.grid.minor = element_blank()) -> beta_plot

  return(beta_plot)

}


## Residuals diagnostics
diagnostics_plot <- function(dt, df_wide, conditioning, associates, phase_name, lead_id) {

  edges <- df_wide %>%
    dplyr::group_by(Block) %>%
    summarise(min = min(Time),
              max = max(Time)) %>%
    ungroup()
  edges_bars <- round(c(edges$min, last(edges$max)), 1)
  edges_line <- geom_vline(xintercept = edges_bars, col = "red", size = 0.75, linetype = "dashed")
  edges_scale <- scale_x_continuous(breaks = edges_bars, labels = edges_bars, expand = c(0.025, 0.025))

  time_min <- round(min(dt$Time))
  time_max <- round(max(dt$Time))
  res_fit <- dt %>%
    dplyr::filter(boot_id == 1)

  res_fit %>%
    ggplot(aes(x = scale(Residuals), y = stat(density))) +
    geom_histogram(position = "identity", fill = "black", alpha = 0.925, bins = 25) +
    scale_x_continuous(breaks = number_ticks(8), expand = c(0, 0)) +
    scale_y_continuous(breaks = number_ticks(8), expand = c(0, 0)) +
    ggtitle(paste0("Channel ", associates, "|", conditioning, " ", phase_name, " Goodness-of-fit"),
            subtitle = paste0("Residuals Histogram / Lag = ", lead_id)) +
    labs(y = "Density", x = "Residual") +
    my_theme + theme(panel.grid.minor = element_blank()) -> res_hist

  res_fit %>% ggplot(aes(sample = scale(Residuals, scale = FALSE))) +
    stat_qq(size = 2.5) +
    stat_qq_line(col = "red") +
    scale_x_continuous(breaks = number_ticks(8), expand = c(0, 0)) +
    scale_y_continuous(breaks = number_ticks(8), expand = c(0, 0)) +
    ggtitle("", subtitle = paste0("Residuals Q-Q plot / Lag = ", lead_id)) +
    labs(y = "Empirical Quantile", x = "Theoretical Quantile") +
    my_theme + theme(panel.grid.minor = element_blank()) -> res_qqplot

  res_fit %>% ggplot() +
    geom_point(aes(x = Time, y = scale(Residuals, scale = FALSE)), size = 2.5) +
    edges_line + edges_scale +
    scale_y_continuous(breaks = number_ticks(8), expand = c(0, .1)) +
    ggtitle("", subtitle = paste0("Residuals on Time / Lag = ", lead_id)) +
    labs(y = "Residuals", x = "Time") +
    my_theme +
    theme(panel.grid.minor = element_blank(),
          plot.margin = unit(c(5.5, 12.5, 5.5, 5.5), "pt"),
          axis.text.x = element_text(size = 12)) -> res_time

  return(cowplot::plot_grid(res_hist, res_qqplot, res_time, ncol = 3))

}
