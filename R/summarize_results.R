## Summarize simulation or other results

###############################################################################
########## Summarize single simulation for all scenarios/models ###############
###############################################################################
simulation_1234_summary <- function() {
  path <- "simulations_1234/"
  model_list <- c("LowFR", "horseshoe", "FIN", "BKMR", "BKMR_group", "CorrQuadReg")
  for (m in model_list) {
    for (scenario in 1:3) {
      # load data for given model and scenario
      curr_path <- paste0(path, m, "/scenario", scenario, "/")
      curr_files <- list.files(path=curr_path)
      curr_path <- paste0(curr_path, curr_files[1])
      curr_data <- readRDS(curr_path)
      
      # summarize the results
      
    }
  }
}

###############################################################################
######## Summarize single simulation for a given scenario and model ###########
###############################################################################
single_sim_summary <- function(model, scenario, print_mixing=FALSE) {
  # load specified data
  path <- "simulations_1234/"
  curr_path <- paste0(path, model, "/scenario", scenario, "/")
  curr_file <- list.files(path=curr_path)
  curr_path <- paste0(curr_path, curr_file[1])
  curr_data <- readRDS(curr_path)
  
  # print summary
  print(paste0("Scenario ", scenario, " ", model, " summary:"))
  print("Accuracy:")
  print(curr_data$accuracy_summary)
  if (print_mixing) {
    print("Mixing:")
    print(curr_data$diagnostic_summary)
  }
  
  return(0)
}








