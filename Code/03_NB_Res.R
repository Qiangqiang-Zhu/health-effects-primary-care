#### Script for computing Moran's I and ACF for negative binomial models

library(tidyverse)
library(data.table)
library(sf)
library(spdep)

# Tidy data ----

load("./Results/Result_NB_Model.RData")  # for the full study period
# load("./Results/Result_NB_Model_v2.RData")

gp_loc <- fread("./Data/GP Practices - Scotland.csv") %>%
  .[, .(PracticeCode = as.character(PracticeCode), Latitude, Longitude)]

result_list <- list()
for (v in ls(pattern = "^nb_")) {
  Medication <- case_when(
      str_extract(v, "bd|ics|rx") == "bd" ~ "Reliever",
      str_extract(v, "bd|ics|rx") == "ics" ~ "Preventer",
      TRUE ~ "Both"
    )

  Model <- paste("Model", str_extract_all(v, "\\d+"))

  result_list[[v]] <- get(v)$results %>% mutate(Medication = Medication, Model = Model)
}

result_df <- bind_rows(result_list) %>%
  mutate(
    Pollutant = factor(
      Pollutant,
      levels = c("NO2", "PM10", "PM25"),
      labels = c(expression(NO[2]), expression(PM[10]), expression(PM[2.5]))
    )
  ) %>%
  left_join(gp_loc, by = "PracticeCode")

# Compute Moran's I ----

moran_calc <- function(sub_data) {
  sub_data <- sub_data %>%
    group_by(Latitude, Longitude) %>%
    summarise(Res = mean(Res), .groups = "drop")

  coords <- data.frame(lon = sub_data$Longitude, lat = sub_data$Latitude)
  nb <- knearneigh(coords, k = 5) %>% knn2nb() %>% make.sym.nb()
  weights <- nb2listw(nb, style = "B")
  moran_test <- moran.mc(sub_data$Res, listw = weights, nsim = 10000)
  moran_stat <- moran_test$statistic
  moran_p1 <- moran_test$p.value
  moran_p2 <- moran.mc(sub_data$Res, listw = weights, nsim = 10000, alternative = "two.sided")$p.value

  data.frame(Moran_I = moran_stat, Moran_P1 = moran_p1, Moran_P2 = moran_p2)
}

moran_df <- result_df %>%
  dplyr::select(-c(Obs, Fit)) %>%
  group_by(Year_Month, Pollutant, Medication, Model) %>%
  nest() %>%
  mutate(results = map(data, moran_calc)) %>%
  unnest(results) %>%
  dplyr::select(-data) %>%
  ungroup()

## Compute ACF(1) ----

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

save(result_df, moran_df, acf_df, file = "./Results/Residual_NB_Model.RData")
# save(result_df, moran_df, acf_df, file = "./Results/Residual_NB_Model_v2.RData")
