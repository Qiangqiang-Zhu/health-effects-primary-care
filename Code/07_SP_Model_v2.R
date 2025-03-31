#### Script for Bayesian spatio-temporal models based on the pre-pandemic

library(tidyverse)
library(data.table)
library(INLA)
library(sf)
library(spdep)
library(doParallel)
library(foreach)

load("./Data/Modelling Data.RData")
data1 <- data %>% filter(Year != 2020)
rm(data)

# Function for half-normal priors in INLA ----

half_normal_tau <- function(sigma){
  return(
    paste("expression:
              sigma = ", sigma, ";
              prec = exp(theta);
              log_dens = -0.5 * log(2 * pi * sigma^2) - 1.5 * theta - 1 / (2 * prec * sigma^2) + theta;
              return(log_dens);",
          sep = '')
  )
}

leroux_prior <- list(
  theta1 = list(prior = half_normal_tau(1.5)),
  theta2 = list(param = c(0, 0.25))
)
ar1_prior <- list(rho = list(param = c(0.25, 0.25)))

# Modelling based on 5NN ----

gp_loc <- filter(gp_loc, PracticeCode %in% unique(data1$PracticeCode))
nb <- st_jitter(st_geometry(gp_loc), factor = 0.0001) %>%
  knearneigh(k = 5) %>%
  knn2nb(sym = TRUE)
nb2INLA("./Data/Scotland.adj", nb)

gp_id <- as.data.frame(gp_loc) %>%
  dplyr::select(PracticeCode) %>%
  mutate(GP = 1:nrow(.))

g <- inla.read.graph(filename = "./Data/Scotland.adj")
data <- left_join(data1, gp_id, by = "PracticeCode")

cl <- makeCluster(16)
registerDoParallel(cl)

grid_5nn <- expand.grid(
  plt = c("no2", "pm10", "pm25"),
  med = c("BD", "ICS", "Rx")
)

model_results <- foreach(
  v = 1:nrow(grid_5nn),
  .export = c("leroux_prior", "ar1_prior", "g"),
  .packages = c("tidyverse", "INLA", "sf"),
  .combine = "c"
) %dopar% {
  i <- grid_5nn$plt[v]
  j <- grid_5nn$med[v]

  formula <- switch(
    j,
    "BD" = as.formula(
      paste0(
        "Number", j, "~ 1 + IncRate + EduAttain + EduPartici + CrimeRate +",
        "HouseOCrat + HouseNCrat + Chinese + Black + Humidity +",
        toupper(i), "+ as.factor(Year) + as.factor(Month) +",
        "f(GP, model = 'besagproper2', graph = g, group = Temp_Trend,",
        "control.group = list(model = 'ar1', hyper = ar1_prior), hyper = leroux_prior)"
      )
    ),
    "ICS" = formula <- as.formula(
      paste0(
        "Number", j, "~ 1 + IncRate + EduAttend + EduAttain + EduPartici +",
        "CrimeRate + HouseNCrat + Chinese + SouthAsian + Black + Humidity + Temperature +",
        toupper(i), "+ as.factor(Year) + as.factor(Month) +",
        "f(GP, model = 'besagproper2', graph = g, group = Temp_Trend,",
        "control.group = list(model = 'ar1', hyper = ar1_prior), hyper = leroux_prior)"
      )
    ),
    as.formula(
      paste0(
        "Number", j, "~ 1 + IncRate + EduAttend + EduAttain + EduPartici +",
        "CrimeRate + HouseNCrat + Chinese + Black + Humidity + Temperature +",
        toupper(i), "+ as.factor(Year) + as.factor(Month) +",
        "f(GP, model = 'besagproper2', graph = g, group = Temp_Trend,",
        "control.group = list(model = 'ar1', hyper = ar1_prior), hyper = leroux_prior)"
      )
    )
  )

  model <- inla(
    formula = formula,
    family = "poisson",
    data = data,
    E = ExpectedNumberRx,
    control.predictor = list(compute = TRUE)
    # control.compute = list(dic = TRUE, waic = TRUE, config = TRUE)
  )

  setNames(list(model), paste0("leroux_", tolower(j), "_", i))
}

for (m in names(model_results)) assign(m, model_results[[m]])

# Save results
save(list = c("data", names(model_results)), file = "./Results/Result_5nn_v2.RData")

stopCluster(cl)

# Modelling based on Mat1 and Mat2 ----

cl <- makeCluster(16)
registerDoParallel(cl)

model_results <- foreach(
  v = c(
    "mat1_eps1", "mat1_eps2", "mat1_eps3", "mat1_eps4", "mat1_eps5",
    "mat2_eps1", "mat2_eps2", "mat2_eps3", "mat2_eps4", "mat2_eps5"
  ),
  .packages = c("tidyverse", "INLA", "sf"),
  .export = c("leroux_prior", "ar1_prior", "data1"),
  .combine = "c"
) %dopar% {
  file_name <- paste0("./Data/Matrices/", v, "_v2.RData")
  load(file_name)

  data <- left_join(data1, adj_id, by = "PracticeCode")
  g <- adj_matrix_sparse

  results <- list()

  for (i in c("no2", "pm10", "pm25")) {
    for (j in c("BD", "ICS", "Rx")) {
      formula <- switch(
        j,
        "BD" = as.formula(
          paste0(
            "Number", j, "~ 1 + IncRate + EduAttain + EduPartici + CrimeRate +",
            "HouseOCrat + HouseNCrat + Chinese + Black + Humidity +",
            toupper(i), "+ as.factor(Year) + as.factor(Month) +",
            "f(GP, model = 'besagproper2', graph = g, group = Temp_Trend,",
            "control.group = list(model = 'ar1', hyper = ar1_prior), hyper = leroux_prior)"
          )
        ),
        "ICS" = formula <- as.formula(
          paste0(
            "Number", j, "~ 1 + IncRate + EduAttend + EduAttain + EduPartici +",
            "CrimeRate + HouseNCrat + Chinese + SouthAsian + Black + Humidity + Temperature +",
            toupper(i), "+ as.factor(Year) + as.factor(Month) +",
            "f(GP, model = 'besagproper2', graph = g, group = Temp_Trend,",
            "control.group = list(model = 'ar1', hyper = ar1_prior), hyper = leroux_prior)"
          )
        ),
        as.formula(
          paste0(
            "Number", j, "~ 1 + IncRate + EduAttend + EduAttain + EduPartici +",
            "CrimeRate + HouseNCrat + Chinese + Black + Humidity + Temperature +",
            toupper(i), "+ as.factor(Year) + as.factor(Month) +",
            "f(GP, model = 'besagproper2', graph = g, group = Temp_Trend,",
            "control.group = list(model = 'ar1', hyper = ar1_prior), hyper = leroux_prior)"
          )
        )
      )

      model <- inla(
        formula = formula,
        family = "poisson",
        data = data,
        E = ExpectedNumberRx,
        control.predictor = list(compute = TRUE)
        # control.compute = list(dic = TRUE, waic = TRUE, config = TRUE)
      )

      results[[paste0("leroux_", tolower(j), "_", i)]] <- model
    }
  }

  for (m in names(results)) assign(m, results[[m]])
  save(list = c("data", names(results)), file = paste0("./Results/Result_", v, "_v2.RData"))
  rm(list = ls())
  gc()
}

stopCluster(cl)
