#### Script for computing Moran's I and ACF for Bayesian models

library(tidyverse)
library(sf)
library(spdep)
library(forecast)
library(data.table)
library(doParallel)
library(foreach)


# Tidy model residuals ----

cl <- makeCluster(16)
registerDoParallel(cl)

result_df <- foreach(
  v = c("5nn", "mat1_eps1", "mat1_eps2", "mat1_eps3", "mat1_eps4", "mat1_eps5",
        "mat2_eps1", "mat2_eps2", "mat2_eps3", "mat2_eps4", "mat2_eps5"),
  .packages = c("tidyverse", "data.table"),
  .combine = rbind
) %dopar% {
  # file_name <- paste0("./Results/Result_", v, ".RData")
  file_name <- paste0("./Results/Result_", v, "_v2.RData")
  load(file_name)

  results <- list()
  index <- 1

  for (i in c("no2", "pm10", "pm25")) {
    for (j in c("BD", "ICS", "Rx")) {
      model_name <- paste0("leroux_", tolower(j), "_", i)
      model <- get(model_name)

      results[[index]] <- data %>%
        dplyr::rename(Obs = !!sym(paste0("Number", j))) %>%
        select(Year, Month, PracticeCode, Obs, ExpectedNumberRx) %>%
        mutate(
          Fit = model$summary.fitted.values[, "0.5quant"] * ExpectedNumberRx,
          Res = (Obs - Fit) / sqrt(Fit),
          Pollutant = toupper(i),
          Medication = j,
          Model = v
        ) %>%
        select(-ExpectedNumberRx)
      index <- index + 1
    }
  }

  bind_rows(results)
}

stopCluster(cl)

result_df <- result_df %>%
  mutate(
    Medication = factor(
      Medication,
      levels = c("BD", "ICS", "Rx"),
      labels = c("Reliever", "Preventer", "Both")
    ),
    Pollutant = factor(
      Pollutant,
      levels = c("NO2", "PM10", "PM25"),
      labels = c(expression(NO[2]), expression(PM[10]), expression(PM[2.5]))
    )
  )

# save(result_df, file = "./Results/Residual_Bayesian_Model.RData")
save(result_df, file = "./Results/Residual_SP_Model_v2.RData")


# Compute Moran's I and ACF(1) ----

load("./Results/Residual_Bayesian_Model_v2.RData")

gp_loc <- fread("./Data/GP Practices - Scotland.csv") %>%
  .[, .(PracticeCode = as.character(PracticeCode), Latitude, Longitude)]

result_df <- result_df %>%
  mutate(Year_Month = sprintf("%04d-%02d", Year, Month)) %>%
  left_join(gp_loc, by = "PracticeCode")

moran_calc <- function(sub_data) {
  sub_data <- sub_data %>%
    group_by(Latitude, Longitude) %>%
    summarise(Res = mean(Res), .groups = "drop")

  coords <- data.frame(lon = sub_data$Longitude, lat = sub_data$Latitude)
  nb <- knearneigh(coords, k = 5) %>% knn2nb() %>% make.sym.nb()
  w(nb, style = "B")
  moran_test <- moran.mc(sub_data$Res, listw = weights, nsim = 10000)
  moran_stat <- moran_test$statistic
  moran_p1 <- moran_test$p.value
  moran_p2 <- moran.mc(sub_data$Res, listw = weights, nsim = 10000, alternative = "two.sided")$p.value

  data.frame(Moran_I = moran_stat, Moran_P1 = moran_p1, Moran_P2 = moran_p2)
}

moran_df <- result_df %>%
  dplyr::select(-c(Obs, Fit)) %>%
  group_by(Year, Month, Pollutant, Medication, Model) %>%
  nest() %>%
  mutate(results = map(data, moran_calc)) %>%
  unnest(results) %>%
  dplyr::select(-data) %>%
  ungroup()

acf_calc <- function(sub_data) {
  acf_result <- acf(sub_data$Res, lag.max = 1, plot = FALSE)
  return(acf_result$acf[2])
}

acf_df <- result_df %>%
  dplyr::select(PracticeCode, Pollutant, Medication, Model, Res) %>%
  group_by(PracticeCode, Pollutant, Medication, Model) %>%
  nest() %>%
  mutate(ACF = map_dbl(data, acf_calc)) %>%
  dplyr::select(-data) %>%
  ungroup()

# save(moran_df, acf_df, file = "./Results/Moran_ACF_Bayesian_Model.RData")
save(moran_df, acf_df, file = "./Results/Moran_ACF_SP_Model_v2.RData")
