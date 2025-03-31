#### Script for visualising results of Bayesian models

library(tidyverse)
library(data.table)
library(ggsci)
library(sf)


# Load and tidy data ----

options(scipen = 10L)
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

load("./Results/Moran_ACF_SP_Model.RData")
moran_df1 <- moran_df

load("./Results/Moran_ACF_SP_Model_v2.RData")
moran_df2 <- moran_df

moran_df <- rbind(
  mutate(moran_df1, Period = "2016-2020"),
  mutate(moran_df2, Period = "2016-2019")
) %>%
  mutate(
    Medication = factor(Medication, levels = c("Reliever", "Preventer", "Both")),
    Period = factor(Period, levels = c("2016-2020", "2016-2019")),
    Matrix = str_extract(Model, "^[^_]+") %>% str_to_title(),
    Threshold = case_when(
      Model == "mat1_eps1" ~ "0.0001",
      Model == "mat1_eps2" ~ "0.0005",
      Model == "mat1_eps3" ~ "0.001",
      Model == "mat1_eps4" ~ "0.005",
      Model == "mat1_eps5" ~ "0.01",
      Model == "mat2_eps1" ~ "0.005",
      Model == "mat2_eps2" ~ "0.01",
      Model == "mat2_eps3" ~ "0.05",
      Model == "mat2_eps4" ~ "0.1",
      Model == "mat2_eps5" ~ "0.5",
      TRUE ~ "-"
    )
  ) %>%
  select(-Model)

# # Tidy the table showing PSAI, GSAI, and MAI ----

moran_tab <- moran_df %>%
  group_by(Pollutant, Medication, Matrix, Period, Threshold) %>%
  summarise(PSAI = sum(Moran_P1 < 0.05) / 60 * 100,
            GSAI = sum(Moran_P2 < 0.05) / 60 * 100,
            MAI = round(mean(abs(Moran_I)), 4),
            .groups = "drop") %>%
  mutate(PSAI = sprintf("%.2f%%", PSAI), GSAI = sprintf("%.2f%%", GSAI)) %>%
  arrange(Pollutant, Medication, Matrix, Threshold, Period)

moran_tab %>% filter(Pollutant == "PM[2.5]") %>% View()

# Moran's I plot for the models based on 5nn ----

moran_5nn_no2 <- moran_df %>%
  filter(Matrix == "5nn", Pollutant == "NO[2]") %>%
  mutate(Year_Month = sprintf("%04d-%02d", Year, Month))


moran_5nn_no2_plot <- moran_5nn_no2 %>%
  ggplot(aes(x = Year_Month, y = Moran_I, colour = Period,
             group = interaction(Period, Medication))) +
  geom_line(linewidth = 0.5) +
  scale_x_discrete(labels = axis_label) +
  labs(title = NULL, x = NULL, y = "Moran's I") +
  scale_color_lancet() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 7),
        legend.margin = margin(t = -5, b = -5, unit = "pt")) +
  geom_text(
    data = moran_5nn_no2 %>%
      group_by(Period, Medication) %>%
      summarise(PSAI = sum(Moran_P1 < 0.05) / 60 * 100, .groups = "drop") %>%
      mutate(vtext = if_else(Period == "2016-2020", 0.04, 0.03)),
    aes(x = 5, y = vtext, label = sprintf("%.2f%%", PSAI), colour = Period),
    vjust = 0, hjust = 0, inherit.aes = FALSE, size = 3, show.legend = FALSE
  ) +
  geom_text(
    data = moran_5nn_no2 %>%
      group_by(Period, Medication) %>%
      summarise(GSAI = sum(Moran_P2 < 0.05) / 60 * 100, .groups = "drop") %>%
      mutate(vtext = if_else(Period == "2016-2020", 0.04, 0.03)),
    aes(x = 20, y = vtext, label = sprintf("%.2f%%", GSAI), colour = Period),
    vjust = 0, hjust = 0, inherit.aes = FALSE, size = 3, show.legend = FALSE
  ) +
  geom_text(
    data = moran_5nn_no2 %>%
      group_by(Period, Medication) %>%
      summarise(MAI = mean(abs(Moran_I)), .groups = "drop") %>%
      mutate(vtext = if_else(Period == "2016-2020", 0.04, 0.03)),
    aes(x = 35, y = vtext, label = sprintf("%.4f", MAI), colour = Period),
    vjust = 0, hjust = 0, inherit.aes = FALSE, size = 3, show.legend = FALSE
  ) +
  annotate("text", x = 5, y = 0.055, label = "PSAI",
           hjust = 0, vjust = 0, size = 4, color = "black") +
  annotate("text", x = 20, y = 0.055, label = "GSAI",
           hjust = 0, vjust = 0, size = 4, color = "black") +
  annotate("text", x = 35, y = 0.055, label = "MAI",
           hjust = 0, vjust = 0, size = 4, color = "black") +
  facet_grid(~ Medication, labeller = label_parsed)
ggsave(moran_5nn_no2_plot, filename = paste0("./Figures/moran_5nn_no2.jpeg"),
       width = 8, height = 4, units = "in", dpi = 300)

moran_5nn_no2 <- moran_df %>%
  filter(Model == "5nn", Pollutant == "NO[2]") %>%
  mutate(Year_Month = sprintf("%04d-%02d", Year, Month))


data_moran_no2 <- data_moran %>% filter(Pollutant == "NO[2]")
data_moran_pm10 <- data_moran %>% filter(Pollutant == "PM[10]")
data_moran_pm25 <- data_moran %>% filter(Pollutant == "PM[2.5]")

# Moran's I plot for the models based on Mat1 ----

moran_mat1_no2 <- moran_df %>%
  filter(Matrix == "Mat1", Pollutant == "NO[2]") %>%
  mutate(Year_Month = sprintf("%04d-%02d", Year, Month))

moran_mat1_no2_plot <- moran_mat1_no2 %>%
  ggplot(aes(x = Year_Month, y = Moran_I, colour = Threshold,
             group = interaction(Period, Threshold, Medication))) +
  geom_line(linewidth = 0.5) +
  scale_x_discrete(labels = axis_label) +
  labs(title = NULL, x = NULL, y = "Moran's I") +
  scale_color_manual(
    values = pal_lancet()(5),
    labels = c(
      expression(epsilon == 0.0001),
      expression(epsilon == 0.0005),
      expression(epsilon == 0.001),
      expression(epsilon == 0.005),
      expression(epsilon == 0.01)
    )
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.margin = margin(t = -5, b = -5, unit = "pt")) +
  geom_text(
    data = moran_mat1_no2 %>%
      group_by(Period, Threshold, Medication) %>%
      summarise(PSAI = sum(Moran_P1 < 0.05) / 60 * 100, .groups = "drop") %>%
      mutate(vtext = case_when(
        Threshold == 0.0001 ~ 0.16,
        Threshold == 0.0005 ~ 0.14,
        Threshold == 0.001 ~ 0.12,
        Threshold == 0.005 ~ 0.10,
        TRUE ~ 0.08
      )),
    aes(x = 5, y = vtext, label = sprintf("%.2f%%", PSAI), colour = Threshold),
    vjust = 0, hjust = 0, inherit.aes = FALSE, size = 3, show.legend = FALSE
  ) +
  geom_text(
    data = moran_mat1_no2 %>%
      group_by(Period, Threshold, Medication) %>%
      summarise(GSAI = sum(Moran_P2 < 0.05) / 60 * 100, .groups = "drop") %>%
      mutate(vtext = case_when(
        Threshold == 0.0001 ~ 0.16,
        Threshold == 0.0005 ~ 0.14,
        Threshold == 0.001 ~ 0.12,
        Threshold == 0.005 ~ 0.10,
        TRUE ~ 0.08
      )),
    aes(x = 15, y = vtext, label = sprintf("%.2f%%", GSAI), colour = Threshold),
    vjust = 0, hjust = 0, inherit.aes = FALSE, size = 3, show.legend = FALSE
  ) +
  geom_text(
    data = moran_mat1_no2 %>%
      group_by(Period, Threshold, Medication) %>%
      summarise(MAI = mean(abs(Moran_I)), .groups = "drop") %>%
      mutate(vtext = case_when(
        Threshold == 0.0001 ~ 0.16,
        Threshold == 0.0005 ~ 0.14,
        Threshold == 0.001 ~ 0.12,
        Threshold == 0.005 ~ 0.10,
        TRUE ~ 0.08
      )),
    aes(x = 25, y = vtext, label = sprintf("%.4f", MAI), colour = Threshold),
    vjust = 0, hjust = 0, inherit.aes = FALSE, size = 3, show.legend = FALSE
  ) +
  annotate("text", x = 5, y = 0.19, label = "PSAI", hjust = 0, vjust = 0, size = 4, color = "black") +
  annotate("text", x = 15, y = 0.19, label = "GSAI", hjust = 0, vjust = 0, size = 4, color = "black") +
  annotate("text", x = 25, y = 0.19, label = "MAI", hjust = 0, vjust = 0, size = 4, color = "black") +
  facet_grid(Medication ~ Period, labeller = label_parsed, scales = "free_x")

ggsave(moran_mat1_no2_plot, filename = paste0("./Figures/moran_mat1_no2.jpeg"),
       width = 8, height = 8, units = "in", dpi = 300)

# Moran's I plot for the models based on Mat2 ----

moran_mat2_no2 <- moran_df %>%
  filter(Matrix == "Mat2", Pollutant == "NO[2]") %>%
  mutate(Year_Month = sprintf("%04d-%02d", Year, Month))

moran_mat2_no2_plot <- moran_mat2_no2 %>%
  ggplot(aes(x = Year_Month, y = Moran_I, colour = Threshold,
             group = interaction(Period, Threshold, Medication))) +
  geom_line(linewidth = 0.5) +
  scale_x_discrete(labels = axis_label) +
  labs(title = NULL, x = NULL, y = "Moran's I") +
  scale_color_manual(
    values = pal_lancet()(5),
    labels = c(
      expression(alpha == 0.005),
      expression(alpha == 0.01),
      expression(alpha == 0.05),
      expression(alpha == 0.1),
      expression(alpha == 0.5)
    )
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.margin = margin(t = -5, b = -5, unit = "pt")) +
  geom_text(
    data = moran_mat2_no2 %>%
      group_by(Period, Threshold, Medication) %>%
      summarise(PSAI = sum(Moran_P1 < 0.05) / 60 * 100, .groups = "drop") %>%
      mutate(vtext = case_when(
        Threshold == 0.005 ~ 0.16,
        Threshold == 0.01 ~ 0.14,
        Threshold == 0.05 ~ 0.12,
        Threshold == 0.1 ~ 0.10,
        TRUE ~ 0.08
      )),
    aes(x = 5, y = vtext, label = sprintf("%.2f%%", PSAI), colour = Threshold),
    vjust = 0, hjust = 0, inherit.aes = FALSE, size = 3, show.legend = FALSE
  ) +
  geom_text(
    data = moran_mat2_no2 %>%
      group_by(Period, Threshold, Medication) %>%
      summarise(GSAI = sum(Moran_P2 < 0.05) / 60 * 100, .groups = "drop") %>%
      mutate(vtext = case_when(
        Threshold == 0.005 ~ 0.16,
        Threshold == 0.01 ~ 0.14,
        Threshold == 0.05 ~ 0.12,
        Threshold == 0.1 ~ 0.10,
        TRUE ~ 0.08
      )),
    aes(x = 15, y = vtext, label = sprintf("%.2f%%", GSAI), colour = Threshold),
    vjust = 0, hjust = 0, inherit.aes = FALSE, size = 3, show.legend = FALSE
  ) +
  geom_text(
    data = moran_mat2_no2 %>%
      group_by(Period, Threshold, Medication) %>%
      summarise(MAI = mean(abs(Moran_I)), .groups = "drop") %>%
      mutate(vtext = case_when(
        Threshold == 0.005 ~ 0.16,
        Threshold == 0.01 ~ 0.14,
        Threshold == 0.05 ~ 0.12,
        Threshold == 0.1 ~ 0.10,
        TRUE ~ 0.08
      )),
    aes(x = 25, y = vtext, label = sprintf("%.4f", MAI), colour = Threshold),
    vjust = 0, hjust = 0, inherit.aes = FALSE, size = 3, show.legend = FALSE
  ) +
  annotate("text", x = 5, y = 0.19, label = "PSAI", hjust = 0, vjust = 0, size = 4, color = "black") +
  annotate("text", x = 15, y = 0.19, label = "GSAI", hjust = 0, vjust = 0, size = 4, color = "black") +
  annotate("text", x = 25, y = 0.19, label = "MAI", hjust = 0, vjust = 0, size = 4, color = "black") +
  facet_grid(Medication ~ Period, labeller = label_parsed, scales = "free_x")
ggsave(moran_mat2_no2_plot, filename = paste0("./Figures/moran_mat2_no2.jpeg"),
       width = 8, height = 8, units = "in", dpi = 300)
