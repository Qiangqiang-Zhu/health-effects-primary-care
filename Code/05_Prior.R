library(tidyverse)
library(INLA)
library(patchwork)

old <- theme_set(
  theme_bw() +
    theme(
      axis.title.x = element_text(size = 14),
      axis.text.x = element_text(size = 12),
      axis.title.y = element_text(size = 14),
      axis.text.y = element_text(size = 12),
      title = element_text(size = 16),
      legend.text = element_text(size = 14),
      legend.title = element_blank(),
      legend.position = "bottom"
    )
)

# Prior for tau/sigma ----

prior_sigma <- function(x, sd = 1) {
  dsigma <- sqrt(2/pi) / sd * exp(-0.5 * x^2 / sd^2)
  return(dsigma)
}

plot_sigma <- function(sd0 = 1) {
  sigma <- seq(0, 6, 0.01)
  pdf_sigma <- prior_sigma(sigma, sd = sd0)
  data.frame(sigma = sigma, pdf = pdf_sigma) %>%
    ggplot(aes(x = sigma, y = pdf)) +
    geom_line(col = "steelblue", linewidth = 0.8) +
    ylim(0, 0.75) +
    labs(title = bquote(paste("PDF for ", sigma)), x = NULL, y = NULL)
}

p1 <- plot_sigma(sd0 = 1.5)

# Prior for rho ----

prior_rho <- function(x, mu = 0, sd = 1) {
  theta <- log(x / (1 - x))
  drho <- 1 / sqrt(2 * pi * sd^2) * exp(-0.5 * ((theta - mu) / sd)^2) / (x * (1 - x))
  return(drho)
}

plot_rho <- function(mean = 0, sd = 1) {
  r_logitrho <- rnorm(n = 1e5, mean = mean, sd = sd)
  r_rho <- exp(r_logitrho) / (1 + exp(r_logitrho))

  rho <- seq(0.01, 0.99, 0.01)
  pdf_rho <- prior_rho(rho, mean, sd)

  ggplot(data.frame(rho = r_rho), aes(x = r_rho)) +
    geom_line(data = data.frame(rho = rho, pdf = pdf_rho), aes(x = rho, y = pdf),
              color = "steelblue", linewidth = 0.8) +
    ylim(0, 2) +
    labs(title = bquote(paste("PDF for ", rho)), x = NULL, y = NULL)
}

p2 <- plot_rho(mean = 0, sd = 2)


# Prior for alpha ----

prior_alpha <- function(x, mu = 0, tau = 1) {
  theta <- log((1 + x) / (1 - x))
  sd <- 1 / sqrt(tau)
  dalpha <- 1 / sqrt(2 * pi * sd^2) * exp(-0.5 * ((theta - mu) / sd)^2) * 2 / (1 - x^2)
  return(dalpha)
}

plot_alpha <- function(mean = 0, sd = 1) {
  r_trans_alpha <- rnorm(n = 1e5, mean = mean, sd = sd)
  r_alpha <- (exp(r_trans_alpha) - 1) / (exp(r_trans_alpha) + 1)

  alpha <- seq(-0.99, 0.99, 0.01)
  pdf_alpha <- prior_alpha(alpha, mean, 1/sd^2)

  ggplot(data.frame(alpha = r_alpha), aes(x = r_alpha)) +
    geom_line(data = data.frame(alpha = alpha, pdf = pdf_alpha), aes(x = alpha, y = pdf),
              color = "steelblue", linewidth = 0.8) +
    ylim(0, 1) +
    labs(title = bquote(paste("PDF for ", alpha)), x = NULL, y = NULL)
}

p3 <- plot_alpha(mean = 0.25, sd = 2)


patch_design <- c("1122
                   #33#")
free(p1) + free(p2) + free(p3) + plot_layout(design = patch_design)
ggsave(filename = "./Figures/priors.jpeg", width = 8, height = 5, units = "in", dpi = 300)
