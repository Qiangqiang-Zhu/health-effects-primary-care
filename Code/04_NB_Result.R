#### Script for tidying and visualising results of negative binomial models

library(tidyverse)
library(data.table)
library(ggsci)
library(sf)
library(spdep)

# Load data ----

load("./Results/Residual_NB_Model.RData")
moran_df1 <- moran_df
acf_df1 <- acf_df

load("./Results/Residual_NB_Model.RData")
moran_df2 <- moran_df
acf_df2 <- acf_df

old <- theme_set(
  theme_bw() +
    theme(
      axis.title.x = element_text(size = 14),
      axis.text.x = element_text(size = 12),
      axis.title.y = element_text(size = 14),
      axis.text.y = element_text(size = 12),
      title = element_text(size = 16),
      strip.text = element_text(size = 12),
      legend.text = element_text(size = 14),
      legend.title = element_blank(),
      legend.position = "bottom"
    )
)

axis_label <- format(
  seq.Date(
    from = as.Date("2016-01-01"),
    to = as.Date("2020-12-01"),
    by = "month"
  ),
  "%Y-%m"
)
for (i in 1:20) {
  axis_label[(i - 1) * 3 + 2] <- ""
  axis_label[(i - 1) * 3 + 3] <- ""
}

## Moran's I plot ##############################################################

data_moran1 <- moran_df1 %>%
  mutate(
    Pollutant = as.factor(Pollutant),
    Medication = factor(
      Medication,
      levels = c("Reliever", "Preventer", "Both")
    ),
    Model = as.factor(
      case_when(Model == "Model 1" ~ "NB 1",
                Model == "Model 2" ~ "NB 2",
                TRUE ~ "NB 3")
    ),
    Period = "2016-2020"
  )

data_moran2 <- moran_df2 %>%
  mutate(
    Pollutant = as.factor(Pollutant),
    Medication = factor(
      Medication,
      levels = c("Reliever", "Preventer", "Both")
    ),
    Model = as.factor(if_else(Model == "Model 1", "NB 1", "NB 3")),
    Period = "2016-2019"
  )

data_moran <- rbind(data_moran1, data_moran2) %>%
  mutate(Period = factor(
    Period,
    levels = c("2016-2020", "2016-2019"),
    labels = c("2016-2020", "2016-2019")
  ))

moran_tab <- data_moran2 %>%
  group_by(Pollutant, Medication, Model) %>%
  summarise(PSAI = sum(Moran_P1 < 0.05) / 60 * 100,
            GSAI = sum(Moran_P2 < 0.05) / 60 * 100,
            MAI = round(mean(abs(Moran_I)), 4),
            .groups = "drop") %>%
  mutate(PSAI = sprintf("%.2f%%", PSAI),
         GSAI = sprintf("%.2f%%", GSAI)) %>%
  arrange(Pollutant, Medication, Model)
View(moran_tab)

data_moran_no2 <- data_moran %>% filter(Pollutant == "NO[2]")

moran_plot <- data_moran_no2 %>%
  ggplot(aes(x = Year_Month, y = Moran_I, colour = Model,
             group = interaction(Period, Model, Medication))) +
  geom_line(linewidth = 0.5) +
  scale_x_discrete(labels = axis_label) +
  labs(title = NULL, x = NULL, y = "Moran's I") +
  scale_color_lancet() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.margin = margin(t = -5, b = -5, unit = "pt")) +
  geom_text(
    data = data_moran_no2 %>%
      group_by(Period, Model, Medication) %>%
      summarise(PSAI = sum(Moran_P1 < 0.05) / 60 * 100, .groups = "drop") %>%
      mutate(vtext = case_when(Model == "NB 1" ~ 0.20, Model == "NB 2" ~ 0.16, TRUE ~ 0.12)),
    aes(x = 5, y = vtext, label = sprintf("%.2f%%", PSAI), colour = Model),
    vjust = 0, hjust = 0,
    inherit.aes = FALSE,
    size = 3,
    show.legend = FALSE
  ) +
  geom_text(
    data = data_moran_no2 %>%
      group_by(Period, Model, Medication) %>%
      summarise(GSAI = sum(Moran_P2 < 0.05) / 60 * 100, .groups = "drop") %>%
      mutate(vtext = case_when(Model == "NB 1" ~ 0.20, Model == "NB 2" ~ 0.16, TRUE ~ 0.12)),
    aes(x = 15, y = vtext, label = sprintf("%.2f%%", GSAI), colour = Model),
    vjust = 0, hjust = 0,
    inherit.aes = FALSE,
    size = 3,
    show.legend = FALSE
  ) +
  geom_text(
    data = data_moran_no2 %>%
      group_by(Period, Model, Medication) %>%
      summarise(MAI = mean(abs(Moran_I)), .groups = "drop") %>%
      mutate(vtext = case_when(Model == "NB 1" ~ 0.20, Model == "NB 2" ~ 0.16, TRUE ~ 0.12)),
    aes(x = 25, y = vtext, label = sprintf("%.4f", MAI), colour = Model),
    vjust = 0, hjust = 0,
    inherit.aes = FALSE,
    size = 3,
    show.legend = FALSE
  ) +
  annotate("text", x = 5, y = 0.25, label = "PSAI",
           hjust = 0, vjust = 0, size = 4, color = "black") +
  annotate("text", x = 15, y = 0.25, label = "GSAI",
           hjust = 0, vjust = 0, size = 4, color = "black") +
  annotate("text", x = 25, y = 0.25, label = "MAI",
           hjust = 0, vjust = 0, size = 4, color = "black") +
  facet_grid(Medication ~ Period, labeller = label_parsed, scales = "free_x")
ggsave(moran_plot, filename = "./Figures/moran_nb_no2.jpeg",
       width = 8, height = 6, units = "in", dpi = 300)

## ACF Plot ####################################################################

data_acf1 <- acf_df1 %>%
  mutate(
    Pollutant = as.factor(Pollutant),
    Medication = factor(Medication, levels = c("Reliever", "Preventer", "Both")),
    Model = as.factor(case_when(Model == "Model 1" ~ "NB 1",
                                Model == "Model 2" ~ "NB 2",
                                TRUE ~ "NB 3"))
  )

acf_plot1 <- data_acf1 %>%
  filter(Pollutant == "NO[2]") %>%
  ggplot(aes(x = PracticeCode, y = ACF, group = interaction(Medication, Model))) +
  geom_point(size = 0.1, colour = "steelblue") +
  scale_x_discrete(labels = NULL) +
  labs(title = NULL, x = NULL, y = "ACF(1)") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "white"),
        legend.margin = margin(t = -5, b = -5, unit = "pt"),
        legend.position = "right") +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  geom_hline(yintercept = 1.96 / sqrt(60), linetype = "dashed", color = "red") +
  geom_hline(yintercept = -1.96 / sqrt(60), linetype = "dashed", color = "red") +
  geom_text(
    data = data_acf1 %>%
      filter(Pollutant == "NO[2]") %>%
      group_by(Model, Medication) %>%
      summarise(sig = sum(abs(ACF) >= 1.96 / sqrt(60)) / n() * 100, .groups = "drop"),
    aes(x = 10, y = -0.8, label = sprintf("%.2f%%", sig)),
    colour = "red",
    vjust = 0, hjust = 0,
    inherit.aes = FALSE,
    size = 3.5,
    show.legend = FALSE
  ) +
  facet_grid(Medication ~ Model)
ggsave(acf_plot1, filename = "./Figures/acf_nb_no2.jpeg",
       width = 8, height = 5, units = "in", dpi = 300)

data_acf2 <- acf_df2 %>%
  mutate(
    Pollutant = as.factor(Pollutant),
    Medication = factor(Medication, levels = c("Reliever", "Preventer", "Both")),
    Model = as.factor(ifelse(Model == "Model 1", "NB 1", "NB 3"))
  )

acf_plot2 <- data_acf2 %>%
  filter(Pollutant == "NO[2]") %>%
  ggplot(aes(x = PracticeCode, y = ACF, group = interaction(Medication, Model))) +
  geom_point(size = 0.1, colour = "steelblue") +
  scale_x_discrete(labels = NULL) +
  labs(title = NULL, x = NULL, y = "ACF(1)") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "white"),
        legend.margin = margin(t = -5, b = -5, unit = "pt"),
        legend.position = "right") +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  geom_hline(yintercept = 1.96 / sqrt(60), linetype = "dashed", color = "red") +
  geom_hline(yintercept = -1.96 / sqrt(60), linetype = "dashed", color = "red") +
  geom_text(
    data = data_acf2 %>%
      filter(Pollutant == "NO[2]") %>%
      group_by(Model, Medication) %>%
      summarise(sig = sum(abs(ACF) >= 1.96 / sqrt(60)) / n() * 100, .groups = "drop"),
    aes(x = 10, y = -0.8, label = sprintf("%.2f%%", sig)),
    colour = "red",
    vjust = 0, hjust = 0,
    inherit.aes = FALSE,
    size = 3.5,
    show.legend = FALSE
  ) +
  facet_grid(Medication ~ Model)
ggsave(acf_plot2, filename = "./Figures/acf_nb_no2_v2.jpeg",
       width = 8, height = 6, units = "in", dpi = 300)
