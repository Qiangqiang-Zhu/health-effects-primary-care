#### Script for exploratory analysis

library(tidyverse)
library(data.table)
library(ggsci)
library(sf)
library(leaflet)
library(leaflegend)
library(patchwork)
library(ggcorrplot)
library(gridExtra)

old <- theme_set(
  theme_bw() +
    theme(
      axis.title.x = element_text(size = 16),
      axis.text.x = element_text(size = 12),
      axis.title.y = element_text(size = 16),
      axis.text.y = element_text(size = 12),
      legend.text = element_text(size = 14),
      legend.position = "bottom"
    )
)

# Prevalence of Asthma and COPD ----

prevalence <- fread("./Data/Average Prevalence and Incidence (PHS).csv")

ggplot(prevalence, aes(x = Age, y = Rate, fill = Sex)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, colour = "black") +
  labs(y = "Rate per 100 Patients") +
  scale_fill_nejm()
ggsave(filename = "./Figures/prevalence_by_sex_age.jpeg",
       width = 8, height = 5, units = "in", dpi = 300)

# Standard prescription ratios ----

geo_hb <- st_read("./Data/SG_NHS_HealthBoards_2019/SG_NHS_HealthBoards_2019.shp") %>%
  st_transform(crs = 4326)

gp_hb <- fread("./Data/GP Practices - Scotland.csv")[, HBName := gsub("^NHS ", "", HBName)]

spr <- fread("./Data/SPR.csv")

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

# Line plot for SPRs over the study period for Scotland.
p1 <- spr %>% group_by(Year, Month) %>%
  summarise(SPR = mean(SPR), .groups = "drop") %>%
  mutate(
    Year = as.factor(Year),
    Month = factor(Month, levels = c(1:12), labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
  ) %>%
  ggplot(aes(x = Month, y = SPR, group = Year, shape = Year, color = Year)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  labs(title = "(a) Lineplot of mean SPR for each month", y = "SPR") +
  scale_y_continuous(labels = scales::comma_format()) +
  theme(title = element_text(size = 10),
        legend.title = element_blank(),
        legend.margin = margin(t = 0)) +
  scale_color_lancet()

geo_spr <- spr[, .(Year, Month, PracticeCode, SPR)] %>%
  group_by(PracticeCode) %>%
  summarise(SPR = mean(SPR), .groups = "drop") %>%
  mutate(SPRGroup = factor(case_when(SPR > 1.1 ~ "> 1.1", SPR < 0.9 ~ "< 0.9", .default = "[0.9,1.1]"),
                           levels = c("> 1.1", "[0.9,1.1]", "< 0.9"))) %>%
  left_join(gp_hb[, .(PracticeCode, Latitude, Longitude)], by = "PracticeCode")

# Spatial map of mean SPRs across Scotland.
spr_colors <- c("< 0.9" = "seagreen", "[0.9,1.1]" = "orange2", "> 1.1" = "red3")
p2 <- ggplot() +
  geom_sf(data = geo_hb, fill = NA, color = "#444444", linewidth = 0.2) +
  geom_point(data = geo_spr, aes(x = Longitude, y = Latitude, color = SPRGroup),
             size = 0.4, alpha = 0.8) +
  scale_color_manual(values = spr_colors, name = "SPR") +
  ggtitle("(b) Map of mean SPR for each surgery") +
  coord_sf(xlim = c(-8.5, -0.5), ylim = c(54.6, 61), expand = FALSE) +
  theme(title = element_text(size = 10),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.spacing.y = unit(0.2, "cm"),
        legend.spacing.x = unit(0, "cm"),
        legend.margin = margin(t = 0),  # tighten legend spacing
        plot.margin = margin(l = 0, r = 5, t = 0, b = 5)) +
  guides(colour = guide_legend(override.aes = list(size = 2)))

p <- p1 + p2 + plot_layout(ncol = 2, widths = c(1.5, 1))
ggsave(p, filename = "./Figures/spr.jpeg", width = 9, height = 5, units = "in", dpi = 300)

# SIMD indicators ----

simd <- st_read("./Data./SIMD 2020/SG_SIMD_2020/SG_SIMD_2020.shp") %>%
  st_transform(crs = 4326)

simd_dt <- as.data.table(simd) %>%
  .[, .(DataZone,
        IncRate = as.numeric(gsub("%", "", IncRate)) * 0.01,
        EmpRate = as.numeric(gsub("%", "", EmpRate)) * 0.01,
        EduAttend = as.numeric(gsub("%", "", EduAttend)) * 0.01, EduAttain,
        EduNoQuals, EduPartici = as.numeric(gsub("%", "", EduPartici)) * 0.01,
        EduUniver = as.numeric(gsub("%", "", EduUniver)) * 0.01,
        GAccPetrol, GAccDTGP, GAccDTPost, GAccDTPsch, GAccDTSsch,
        GAccDTRet, GAccPTGP, GAccPTPost, GAccPTRet,
        GAccBrdbnd = as.numeric(gsub("%", "", GAccBrdbnd)) * 0.01,
        CrimeRate = CrimeRate / 10000,
        HouseOCrat = as.numeric(gsub("%", "", HouseOCrat)) * 0.01,
        HouseNCrat = as.numeric(gsub("%", "", HouseNCrat)) * 0.01)]

# Correlation plot of the SIMD indicators.
simd_dt[, -1] %>%
  select(!starts_with("GAcc")) %>%
  cor() %>%
  ggcorrplot(method = "square", type = "lower", lab = TRUE, tl.srt = 45,
             lab_size = 4, colors = c("#6D9EC1", "white", "#E46726"),
             show.legend = TRUE)
ggsave(filename = "./Figures/simd_corr.jpeg",
       width = 8, height = 5, units = "in", dpi = 300)

# Density plots of the SIMD indicators.
hist_plot <- function(data, variable) {
  ggplot(data, aes(x = .data[[variable]])) +
    geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "white", color = "black") +
    geom_density(color = "steelblue", lwd = 0.5, fill = "steelblue", alpha = 0.5) +
    theme_bw()
}
variables <- c("IncRate", "EduAttend", "EduAttain", "EduPartici",
               "EduUniver", "CrimeRate", "HouseOCrat", "HouseNCrat")
simd_density <- grid.arrange(grobs = lapply(variables, function(var) hist_plot(simd_dt, var)),
                             ncol = 3, nrow = 3)
ggsave(simd_density, filename = "./Figures/simd_density.jpeg",
       width = 8, height = 6, units = "in", dpi = 300)

summary(simd_dt[EduAttend > 0, EduAttend])
summary(simd_dt[EduAttain > 0, EduAttain])

simd_dt <- simd_dt %>%
  select(DataZone, IncRate, EduAttend, EduAttain, EduPartici, EduUniver,
         CrimeRate, HouseOCrat, HouseNCrat) %>%
  mutate(across(c(EduAttend, EduAttain), ~ if_else(. == 0, min(.[. > 0]), .)))

fwrite(simd_dt, "./Data/SIMD.csv")
