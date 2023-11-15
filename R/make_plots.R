## Functions to make plots for paper
library(tidyverse)

###############################################################################
############ Plot 95% intervals and true values for main effects ##############
###############################################################################
main_effect_sim_plot <- function(model, scenario, y_lim=NULL) {
  # get true values of coefficients
  if (scenario==1) {
    true_coefs <- simulate_scenario1(random_seed=1234)
  }
  else if (scenario==2) {
    true_coefs <- simulate_scenario2(random_seed=1234)
  }
  else if (scenario==3) {
    true_coefs <- simulate_scenario3(random_seed=1234)
  }
  else {
    print("Must specify data generation scenario as one of 1, 2, or 3.")
    return(NULL)
  }
  main_effect_true <- true_coefs$alpha
  
  # load fitted model
  path <- paste0("simulations_1234/", model, "/scenario", scenario, "/",
                 model, "_scenario", scenario, "_seed1234.rds")
  data <- readRDS(path)
  
  # get intervals and true values
  main_effect_lower <- c()
  main_effect_upper <- c()
  for (i in 1:30) {
    main_effect_lower <- c(main_effect_lower, quantile(data$post_samples$alpha[,i], probs=c(0.025)))
    main_effect_upper <- c(main_effect_upper, quantile(data$post_samples$alpha[,i], probs=c(0.975)))
  }
  
  # create df
  labels <- c("x_11", "x_12", "x_13",
              "x_21", "x_22", "x_23",
              "x_31", "x_32", "x_33",
              "x_41", "x_42", "x_43",
              "x_51", "x_52", "x_53",
              "x_61", "x_62", "x_63",
              "x_71", "x_72", "x_73",
              "x_81", "x_82", "x_83",
              "x_91", "x_92", "x_93",
              "x_10,1", "x_10,2", "x_10,3")
  main_effect_interval_df <- data.frame(labels, main_effect_lower, main_effect_upper, main_effect_true)
  main_effect_interval_df$labels <- factor(main_effect_interval_df$labels,
                                           levels=rev(main_effect_interval_df$labels))
  
  # capitalize "horseshoe"
  if (model == "horseshoe") {
    model <- "Horseshoe"
  }
  
  # make plot
  main_plot <- ggplot(data=main_effect_interval_df, aes(x=labels, y=main_effect_true)) +
    geom_errorbar(aes(ymin=main_effect_lower, ymax=main_effect_upper)) +
    geom_hline(yintercept=0, color="black") +
    geom_point(color="red") +
    coord_flip() +
    labs(title=paste0("Main effects (", model, ", Scenario ", scenario, ")"),
         y="", x="") +
    theme(plot.title = element_text(size=20),
          axis.text.x = element_text(size=13),
          axis.text.y = element_text(size=10))
  
  if (!is.null(y_lim)) {
    main_plot <- main_plot + ylim(y_lim[1], y_lim[2])
  }
  
  return(main_plot)
}


###############################################################################
####################### Make main effect trace plot ###########################
###############################################################################
main_effect_traceplot <- function(model, scenario, Exp, t, y_lim=NULL, cex.main=4.5, cex.axis=2, cex.lab=1.4) {
  # get true values of coefficients
  if (scenario==1) {
    true_coefs <- simulate_scenario1(random_seed=1234)
  }
  else if (scenario==2) {
    true_coefs <- simulate_scenario2(random_seed=1234)
  }
  else if (scenario==3) {
    true_coefs <- simulate_scenario3(random_seed=1234)
  }
  else {
    print("Must specify data generation scenario as one of 1, 2, or 3.")
    return(NULL)
  }
  main_effect_true <- true_coefs$alpha[(Exp-1)*3 + t]
  
  # load fitted model
  path <- paste0("simulations_1234/", model, "/scenario", scenario, "/",
                 model, "_scenario", scenario, "_seed1234.rds")
  data <- readRDS(path)
  alpha_samps <- data$post_samples$alpha[,(Exp-1)*3 + t]
  
  # capitalize "horseshoe"
  if (model == "horseshoe") {
    model <- "Horseshoe"
  }
  
  # make plot
  #title = paste0("Samples of x_", Exp, ",", t, " (", model, ", Scenario ", scenario, ")")
  title = model
  y_label = paste0("alpha_", Exp, ",", t)
  if (!is.null(y_lim)) {
    plot(alpha_samps, type="l", ylim=y_lim, main=title, ylab=y_label,
         cex.main=cex.main, cex.axis=cex.axis, cex.lab=cex.lab)
  }
  else {
    plot(alpha_samps, type="l", main=title)
  }
  abline(h=main_effect_true, col="red", lwd=4)
}


###############################################################################
######### Plot 95% intervals and true values for cumulative effects ###########
###############################################################################
cum_effect_sim_plot <- function(model, scenario, y_lim=NULL) {
  # get true values of coefficients
  if (scenario==1) {
    true_coefs <- simulate_scenario1(random_seed=1234)
  }
  else if (scenario==2) {
    true_coefs <- simulate_scenario2(random_seed=1234)
  }
  else if (scenario==3) {
    true_coefs <- simulate_scenario3(random_seed=1234)
  }
  else {
    print("Must specify data generation scenario as one of 1, 2, or 3.")
    return(NULL)
  }
  main_effect_true <- true_coefs$alpha
  
  cum_effect_true <- main_effect_true[seq(1,30,3)] +
    main_effect_true[seq(2,30,3)] +
    main_effect_true[seq(3,30,3)]
  
  # load fitted model
  path <- paste0("simulations_1234/", model, "/scenario", scenario, "/",
                 model, "_scenario", scenario, "_seed1234.rds")
  data <- readRDS(path)
  
  # get intervals and true values
  main_effect_samples <- data$post_samples$alpha
  cum_effect_samples <- main_effect_samples[,seq(1,30,3)] +
    main_effect_samples[,seq(2,30,3)] +
    main_effect_samples[,seq(3,30,3)]
  
  cum_effect_lower <- c()
  cum_effect_upper <- c()
  for (i in 1:10) {
    cum_effect_lower <- c(cum_effect_lower, quantile(cum_effect_samples[,i], probs=c(0.025)))
    cum_effect_upper <- c(cum_effect_upper, quantile(cum_effect_samples[,i], probs=c(0.975)))
  }
  
  # create df
  labels <- c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10")
  cum_effect_interval_df <- data.frame(labels, cum_effect_lower, cum_effect_upper, cum_effect_true)
  cum_effect_interval_df$labels <- factor(cum_effect_interval_df$labels,
                                           levels=rev(cum_effect_interval_df$labels))
  
  # capitalize "horseshoe"
  if (model == "horseshoe") {
    model <- "Horseshoe"
  }
  
  # make plot
  main_plot <- ggplot(data=cum_effect_interval_df, aes(x=labels, y=cum_effect_true)) +
    geom_errorbar(aes(ymin=cum_effect_lower, ymax=cum_effect_upper)) +
    geom_hline(yintercept=0, color="black") +
    geom_point(color="red") +
    coord_flip() +
    labs(title=paste0("Cumulative effects (", model, ", Scenario ", scenario, ")"),
         y="", x="") +
    theme(plot.title = element_text(size=20),
          axis.text.x = element_text(size=13),
          axis.text.y = element_text(size=10))
  
  if (!is.null(y_lim)) {
    main_plot <- main_plot + ylim(y_lim[1], y_lim[2])
  }
  
  return(main_plot)
}



