#### Script for the exploratory analysis using negative binomial models

library(tidyverse)
library(data.table)
library(sf)
library(MASS)
library(mgcv)
library(gratia)
library(broom)
library(spdep)

# Build negative binomial models based on the full period 2016-2020 ----

load("./Data/Modelling Data.RData")
data <- mutate(data, PracticeCode = as.factor(PracticeCode))

air_pollutants <- c("NO2", "PM10", "PM25")

nb_func <- function(formula) {
  estimates <- list()
  results <- list()

  for (pollutant in air_pollutants) {
    model_formula <- update(formula, as.formula(paste(". ~ . +", pollutant)))

    model <- bam(
      formula = model_formula,
      data = data,
      offset = log(ExpectedNumberRx),
      family = nb(link = "log"),
      method = "fREML",
      discrete = TRUE
    )

    estimates[[pollutant]] <- tidy(model, exponentiate = TRUE, conf.int = TRUE,
                                   parametric = TRUE) %>%
      filter(term == pollutant) %>%
      mutate(across(where(is.numeric), ~ round(., 4))) %>%
      rename(Pollutant = term, RR = estimate) %>%
      mutate(`95% CI` = paste0("(", conf.low, ", ", conf.high, ")"),
             AIC = AIC(model),
             BIC = BIC(model)) %>%
      dplyr::select(Pollutant, RR, `95% CI`, AIC, BIC)

    results[[pollutant]] <- data %>%
      rename(Obs = NumberRx) %>%
      dplyr::select(Year_Month, PracticeCode, Obs) %>%
      mutate(Fit = fitted(model), Res = residuals(model, type = "pearson"),
             Pollutant = pollutant)
  }

  estimates <- bind_rows(estimates)
  results <- bind_rows(results)

  return(list(estimates = estimates, results = results))
}

## Stepwise selection for all medication model ----

### Select confounders for both medications ----

# Initial model
nb_rx <- bam(
  NumberRx ~ IncRate + EduAttend + EduAttain + EduPartici + EduUniver +
    CrimeRate + HouseOCrat + HouseNCrat + Chinese + SouthAsian + Black +
    Humidity + Temperature +
    s(PracticeCode, bs = "re"),
  data = data,
  offset = log(ExpectedNumberRx),
  family = nb(link = "log"),
  method = "fREML",
  discrete = TRUE
)
summary(nb_rx)

# Model with the selected covariates
nb_rx <- bam(
  NumberRx ~ IncRate + EduAttend + EduAttain + EduPartici + CrimeRate +
    HouseNCrat + Chinese + Black + Humidity + Temperature +
    s(PracticeCode, bs = "re"),
  data = data,
  offset = log(ExpectedNumberRx),
  family = nb(link = "log"),
  method = "fREML",
  discrete = TRUE
)
summary(nb_rx)

### Add single pollutant and the selected confounders ----

# NB 1
nb_rx1 <- nb_func(
  NumberRx ~ IncRate + EduAttend + EduAttain + EduPartici + CrimeRate +
    HouseNCrat + Chinese + Black + Humidity + Temperature +
    as.factor(Year) + as.factor(Month) +
    s(PracticeCode, bs = "re")
)

# NB 2
nb_rx2 <- nb_func(
  NumberRx ~ IncRate + EduAttend + EduAttain + EduPartici + CrimeRate +
    HouseNCrat + Chinese + Black + Humidity + Temperature +
    as.factor(Year) + as.factor(Month) +
    as.factor(Temp_Trend == 51) + as.factor(Temp_Trend == 52) +
    as.factor(Temp_Trend == 53) + as.factor(Temp_Trend == 54) +
    s(PracticeCode, bs = "re")
)

# NB 3
nb_rx3 <- nb_func(
  NumberRx ~ IncRate + EduAttend + EduAttain + EduPartici + CrimeRate +
    HouseNCrat + Chinese + Black + Humidity + Temperature +
    as.factor(Temp_Trend) +
    s(PracticeCode, bs = "re")
)

## Stepwise selection for BD model ----

## Select confounders for relievers

# Initial model
nb_bd <- bam(
  NumberBD ~ IncRate + EduAttend + EduAttain + EduPartici + EduUniver +
    CrimeRate + HouseOCrat + HouseNCrat + Chinese + SouthAsian + Black +
    Humidity + Temperature +
    s(PracticeCode, bs = "re"),
  data = data,
  offset = log(ExpectedNumberRx),
  family = nb(link = "log"),
  method = "fREML",
  discrete = TRUE
)
summary(nb_bd)

# Model with the selected covariates
nb_bd <- bam(
  NumberBD ~ IncRate + EduAttain + EduPartici + CrimeRate + HouseOCrat +
    HouseNCrat + Chinese + Black + Humidity + Temperature +
    s(PracticeCode, bs = "re"),
  data = data,
  offset = log(ExpectedNumberRx),
  family = nb(link = "log"),
  method = "fREML",
  discrete = TRUE
)
summary(nb_bd)

### Add single pollutant and selected confounders ----

# NB 1
nb_bd1 <- nb_func(
  NumberBD ~ IncRate + EduAttain + EduPartici + CrimeRate + HouseOCrat +
    HouseNCrat + Chinese + Black + Humidity + Temperature +
    as.factor(Year) + as.factor(Month) +
    s(PracticeCode, bs = "re")
)

# NB 2
nb_bd2 <- nb_func(
  NumberBD ~ IncRate + EduAttain + EduPartici + CrimeRate + HouseOCrat +
    HouseNCrat + Chinese + Black + Humidity + Temperature +
    as.factor(Year) + as.factor(Month) +
    as.factor(Temp_Trend == 51) + as.factor(Temp_Trend == 52) +
    as.factor(Temp_Trend == 53) + as.factor(Temp_Trend == 54) +
    s(PracticeCode, bs = "re")
)

# NB 3
nb_bd3 <- nb_func(
  NumberBD ~ IncRate + EduAttain + EduPartici + CrimeRate + HouseOCrat +
    HouseNCrat + Chinese + Black + Humidity + Temperature +
    as.factor(Temp_Trend) +
    s(PracticeCode, bs = "re")
)

## Stepwise selection for ICS model ----

### Select confounders for preventers ----

# Initial model
nb_ics <- bam(
  NumberICS ~ IncRate + EduAttend + EduAttain + EduPartici + EduUniver +
    CrimeRate + HouseOCrat + HouseNCrat + Chinese + SouthAsian + Black +
    Humidity + Temperature +
    s(PracticeCode, bs = "re"),
  data = data,
  offset = log(ExpectedNumberRx),
  family = nb(link = "log"),
  method = "fREML",
  discrete = TRUE
)
summary(nb_ics)

# Model with selected covariates
nb_ics <- bam(
  NumberICS ~ IncRate + EduAttend + EduAttain + EduPartici + CrimeRate +
    HouseOCrat + HouseNCrat + Chinese + Black + Humidity + Temperature +
    s(PracticeCode, bs = "re"),
  data = data,
  offset = log(ExpectedNumberRx),
  family = nb(link = "log"),
  method = "fREML",
  discrete = TRUE
)
summary(nb_ics)

### Add single pollutant and selected confounders ----

# NB 1
nb_ics1 <- nb_func(
  NumberICS ~ IncRate + EduAttend + EduAttain + EduPartici + CrimeRate +
    HouseOCrat + HouseNCrat + Chinese + Black + Humidity + Temperature +
    as.factor(Year) + as.factor(Month) +
    s(PracticeCode, bs = "re")
)

# NB 2
nb_ics2 <- nb_func(
  NumberICS ~ IncRate + EduAttend + EduAttain + EduPartici + CrimeRate +
    HouseOCrat + HouseNCrat + Chinese + Black + Humidity + Temperature +
    as.factor(Year) + as.factor(Month) +
    as.factor(Temp_Trend == 51) + as.factor(Temp_Trend == 52) +
    as.factor(Temp_Trend == 53) + as.factor(Temp_Trend == 54) +
    s(PracticeCode, bs = "re")
)

# NB 3
nb_ics3 <- nb_func(
  NumberICS ~ IncRate + EduAttend + EduAttain + EduPartici + CrimeRate +
    HouseOCrat + HouseNCrat + Chinese + Black + Humidity + Temperature +
    as.factor(Temp_Trend) +
    s(PracticeCode, bs = "re")
)

save(nb_rx1, nb_rx2, nb_rx3, nb_bd1, nb_bd2, nb_bd3,
     nb_ics1, nb_ics2, nb_ics3,
     file = "./Results/Result_NB_Model.RData")


# Build negative binomial models based on the prependemic 2016-2010 ----

rm(list = ls())

load("./Data/Modelling Data.RData")
data <- mutate(data, PracticeCode = as.factor(PracticeCode)) %>%
  filter(Year != 2020)

air_pollutants <- c("NO2", "PM10", "PM25")

nb_func <- function(formula) {
  estimates <- list()
  results <- list()

  for (pollutant in air_pollutants) {
    model_formula <- update(formula, as.formula(paste(". ~ . +", pollutant)))

    model <- bam(
      formula = model_formula,
      data = data,
      offset = log(ExpectedNumberRx),
      family = nb(link = "log"),
      method = "fREML",
      discrete = TRUE
    )

    estimates[[pollutant]] <- tidy(model, exponentiate = TRUE, conf.int = TRUE,
                                   parametric = TRUE) %>%
      filter(term == pollutant) %>%
      mutate(across(where(is.numeric), ~ round(., 4))) %>%
      rename(Pollutant = term, RR = estimate) %>%
      mutate(`95% CI` = paste0("(", conf.low, ", ", conf.high, ")"),
             AIC = AIC(model),
             BIC = BIC(model)) %>%
      dplyr::select(Pollutant, RR, `95% CI`, AIC, BIC)

    results[[pollutant]] <- data %>%
      rename(Obs = NumberRx) %>%
      dplyr::select(Year_Month, PracticeCode, Obs) %>%
      mutate(Fit = fitted(model), Res = residuals(model, type = "pearson"),
             Pollutant = pollutant)
  }

  estimates <- bind_rows(estimates)
  results <- bind_rows(results)

  return(list(estimates = estimates, results = results))
}

## Stepwise selection for all medication model ----

### Select confounders for both medications ----

# Initial model
nb_rx <- bam(
  NumberRx ~ IncRate + EduAttend + EduAttain + EduPartici + EduUniver +
    CrimeRate + HouseOCrat + HouseNCrat + Chinese + SouthAsian + Black +
    Humidity + Temperature +
    s(PracticeCode, bs = "re"),
  data = data,
  offset = log(ExpectedNumberRx),
  family = nb(link = "log"),
  method = "fREML",
  discrete = TRUE
)
summary(nb_rx)

# Model with the selected covariates
nb_rx <- bam(
  NumberRx ~ IncRate + EduAttend + EduAttain + EduPartici + CrimeRate +
    HouseNCrat + Chinese + Black + Humidity + Temperature +
    s(PracticeCode, bs = "re"),
  data = data,
  offset = log(ExpectedNumberRx),
  family = nb(link = "log"),
  method = "fREML",
  discrete = TRUE
)
summary(nb_rx)

### Add single pollutant and the selected confounders ----

# NB 1
nb_rx1 <- nb_func(
  NumberRx ~ IncRate + EduAttend + EduAttain + EduPartici + CrimeRate +
    HouseNCrat + Chinese + Black + Humidity + Temperature +
    as.factor(Year) + as.factor(Month) +
    s(PracticeCode, bs = "re")
)

# NB 3
nb_rx3 <- nb_func(
  NumberRx ~ IncRate + EduAttend + EduAttain + EduPartici + CrimeRate +
    HouseNCrat + Chinese + Black + Humidity + Temperature +
    as.factor(Temp_Trend) +
    s(PracticeCode, bs = "re")
)

## Stepwise selection for BD model ----

## Select confounders for relievers

# Initial model
nb_bd <- bam(
  NumberBD ~ IncRate + EduAttend + EduAttain + EduPartici + EduUniver +
    CrimeRate + HouseOCrat + HouseNCrat + Chinese + SouthAsian + Black +
    Humidity + Temperature +
    s(PracticeCode, bs = "re"),
  data = data,
  offset = log(ExpectedNumberRx),
  family = nb(link = "log"),
  method = "fREML",
  discrete = TRUE
)
summary(nb_bd)

# Model with the selected covariates
nb_bd <- bam(
  NumberBD ~ IncRate + EduAttain + EduPartici + CrimeRate + HouseOCrat +
    HouseNCrat + Chinese + Black + Humidity +
    s(PracticeCode, bs = "re"),
  data = data,
  offset = log(ExpectedNumberRx),
  family = nb(link = "log"),
  method = "fREML",
  discrete = TRUE
)
summary(nb_bd)

### Add single pollutant and selected confounders ----

# NB 1
nb_bd1 <- nb_func(
  NumberBD ~ IncRate + EduAttain + EduPartici + CrimeRate + HouseOCrat +
    HouseNCrat + Chinese + Black + Humidity +
    as.factor(Year) + as.factor(Month) +
    s(PracticeCode, bs = "re")
)

# NB 3
nb_bd3 <- nb_func(
  NumberBD ~ IncRate + EduAttain + EduPartici + CrimeRate + HouseOCrat +
    HouseNCrat + Chinese + Black + Humidity +
    as.factor(Temp_Trend) +
    s(PracticeCode, bs = "re")
)

## Stepwise selection for ICS model ----

### Select confounders for preventers ----

# Initial model
nb_ics <- bam(
  NumberICS ~ IncRate + EduAttend + EduAttain + EduPartici + EduUniver +
    CrimeRate + HouseOCrat + HouseNCrat + Chinese + SouthAsian + Black +
    Humidity + Temperature +
    s(PracticeCode, bs = "re"),
  data = data,
  offset = log(ExpectedNumberRx),
  family = nb(link = "log"),
  method = "fREML",
  discrete = TRUE
)
summary(nb_ics)

# Model with selected covariates
nb_ics <- bam(
  NumberICS ~ IncRate + EduAttend + EduAttain + EduPartici + CrimeRate +
    HouseNCrat + Chinese + SouthAsian + Black + Humidity + Temperature +
    s(PracticeCode, bs = "re"),
  data = data,
  offset = log(ExpectedNumberRx),
  family = nb(link = "log"),
  method = "fREML",
  discrete = TRUE
)
summary(nb_ics)

### Add single pollutant and selected confounders ----

# NB 1
nb_ics1 <- nb_func(
  NumberICS ~ IncRate + EduAttend + EduAttain + EduPartici + CrimeRate +
    HouseNCrat + Chinese + SouthAsian + Black + Humidity + Temperature +
    as.factor(Year) + as.factor(Month) +
    s(PracticeCode, bs = "re")
)

# NB 3
nb_ics3 <- nb_func(
  NumberICS ~ IncRate + EduAttend + EduAttain + EduPartici + CrimeRate +
    HouseNCrat + Chinese + SouthAsian + Black + Humidity + Temperature +
    as.factor(Temp_Trend) +
    s(PracticeCode, bs = "re")
)

save(nb_rx1, nb_rx3, nb_bd1, nb_bd3, nb_ics1, nb_ics3,
     file = "./Results/Result_NB_Model_v2.RData")
