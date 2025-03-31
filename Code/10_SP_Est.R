library(tidyverse)
library(INLA)
library(data.table)
library(doParallel)
library(foreach)
library(sf)

cl <- makeCluster(16)
registerDoParallel(cl)

est_mat <- foreach(
  v = c(
    "5nn", "5nn_v2",
    "mat1_eps1", "mat1_eps2", "mat1_eps3", "mat1_eps4", "mat1_eps5",
    "mat2_eps1", "mat2_eps2", "mat2_eps3", "mat2_eps4", "mat2_eps5",
    "mat1_eps1_v2", "mat1_eps2_v2", "mat1_eps3_v2", "mat1_eps4_v2", "mat1_eps5_v2",
    "mat2_eps1_v2", "mat2_eps2_v2", "mat2_eps3_v2", "mat2_eps4_v2", "mat2_eps5_v2"
  ),
  .packages = c("tidyverse", "INLA", "sf"),
  .combine = rbind
) %dopar% {
  load(paste0("./Results/Result_", v, ".RData"))

  results <- list()
  index <- 1
  for (i in c("no2", "pm10", "pm25")) {
    for (j in c("BD", "ICS", "Rx")) {
      model <- get(paste0("leroux_", tolower(j), "_", i))
      beta <- model$marginals.fixed[[toupper(i)]]

      results[[index]] <- model$summary.fixed[toupper(i), c("0.5quant", "0.025quant", "0.975quant")] %>%
        exp() %>%
        rownames_to_column(var = "Pollutant") %>%
        mutate(
          Medication = switch(j, "BD" = "Reliever", "ICS" = "Preventer", "Both"),
          Mat = v,
          EP = 1 - inla.pmarginal(0, beta)
        )
      index <- index + 1
    }
  }
  rm(leroux_bd_no2, leroux_ics_no2, leroux_rx_no2,
     leroux_bd_pm10, leroux_ics_pm10, leroux_rx_pm10,
     leroux_bd_pm25, leroux_ics_pm25, leroux_rx_pm25,
     adj_id, adj_matrix_sparse)
  do.call("rbind", results)
}

save(est_mat, file = "./Results/Result_SP_Est.RData")

stopCluster(cl)

################################################################################

load("./Results/Result_SP_Est.RData")

load("./Data/Modelling Data.RData")
pollutant_sd <- c("NO2" = sd(data$NO2), "PM10" = sd(data$PM10), "PM25" = sd(data$PM25))

est_tab <- est_mat %>%
  mutate(Medication = factor(Medication, levels = c("Reliever", "Preventer", "Both"))) %>%
  mutate(
    Matrix = str_extract(Mat, "^[^_]+") %>% str_to_title(),
    Period = if_else(str_detect(Mat, "v2"), "2016-2019", "2016-2020"),
    M = str_remove(Mat, "_v2"),
    Threshold = case_when(
      M == "mat1_eps1" ~ "0.0001",
      M == "mat1_eps2" ~ "0.0005",
      M == "mat1_eps3" ~ "0.001",
      M == "mat1_eps4" ~ "0.005",
      M == "mat1_eps5" ~ "0.01",
      M == "mat2_eps1" ~ "0.005",
      M == "mat2_eps2" ~ "0.01",
      M == "mat2_eps3" ~ "0.05",
      M == "mat2_eps4" ~ "0.1",
      M == "mat2_eps5" ~ "0.5",
      TRUE ~ "-"
    ),
    SD = pollutant_sd[Pollutant],
    RR = round(`0.5quant`^SD, 4),
    EP = round(EP, 4),
    CI = paste0("(", round(`0.025quant`^SD, 4), ", ", round(`0.975quant`^SD, 4), ")")
  ) %>%
  arrange(Pollutant, Medication, Matrix, Threshold) %>%
  select(Pollutant, Medication, Matrix, Period, Threshold, RR, EP, CI)

est_tab %>% filter(Pollutant == "NO2") %>% View()

