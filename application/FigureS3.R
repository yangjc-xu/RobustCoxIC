## ---------------------------------------------------------------------------
## Figure S3: residual processes for the endpoint of hypertension.
##
## Run from the root directory of the reproducibility materials, after
## application/gof_test.R has been run for the eight tests shown below:
##   Rscript application/FigureS3.R
##
## Input:  application/output/ff_{pc1,pc4,pc7,pc9}_hypertension.RData,
##         application/output/prop_{sex,england,pc4,pc9}_hypertension.RData
## Output: application/output/FigureS3.pdf and application/output/FigureS3.eps
## ---------------------------------------------------------------------------

if (!file.exists("application/functions.R")) {
  stop("Please run this script from the root directory of the reproducibility materials.")
}
library(tidyverse)
library(ggpubr)
source("application/functions.R")

## Functional forms of covariates
p_hypertension_pc1 = gof_panel(
  "hypertension", "ff", "pc1",
  title = "Hypertension: functional form for PC1", xlab = "PC1",
  xlim = c(-17,-5), xbreaks = seq(-17,-5,2),
  ylim = c(-1,1.5), ybreaks = c(seq(-1,0,0.5),seq(0.5,1.5,0.5)),
  text_x = -10, text_y = -0.75, aspect = 4/8, sim_rows = 201:220)

p_hypertension_pc4 = gof_panel(
  "hypertension", "ff", "pc4",
  title = "Hypertension: functional form for PC4", xlab = "PC4",
  xlim = c(-10,14), xbreaks = seq(-10,14,3),
  ylim = c(-2,0.8), ybreaks = c(seq(-2,0,0.4),seq(0.4,0.8,0.4)),
  text_x = 5, text_y = -1.7, aspect = 4/8)

p_hypertension_pc7 = gof_panel(
  "hypertension", "ff", "pc7",
  title = "Hypertension: functional form for PC7", xlab = "PC7",
  xlim = c(-10,10), xbreaks = seq(-10,10,2),
  ylim = c(-1,0.6), ybreaks = c(seq(-1,0,0.2),seq(0.2,0.6,0.2)),
  text_x = 2.5, text_y = -0.8, aspect = 4/8)

p_hypertension_pc9 = gof_panel(
  "hypertension", "ff", "pc9",
  title = "Hypertension: functional form for PC9", xlab = "PC9",
  xlim = c(-15,9), xbreaks = seq(-15,9,3),
  ylim = c(-2.5,1), ybreaks = c(seq(-2.5,0,0.5),seq(0.5,1,0.5)),
  text_x = -12, text_y = -2, aspect = 4/8)

## Proportional hazards assumption
p_hypertension_prop_sex = gof_panel(
  "hypertension", "prop", "sex",
  title = "Hypertension: proportionality for sex", xlab = "Age at onset",
  xlim = c(38,85), xbreaks = seq(40,85,5),
  ylim = c(-1,3), ybreaks = c(seq(-1,0,0.5),seq(0.5,3,0.5)),
  text_x = 65, text_y = -0.6, aspect = 4/8)

p_hypertension_prop_england = gof_panel(
  "hypertension", "prop", "england",
  title = "Hypertension: proportionality for birthplace", xlab = "Age at onset",
  xlim = c(38,85), xbreaks = seq(40,85,5),
  ylim = c(-0.5,0.3), ybreaks = c(seq(-0.5,0,0.1),seq(0.1,0.3,0.1)),
  text_x = 65, text_y = -0.4, aspect = 4/8)

p_hypertension_prop_pc4 = gof_panel(
  "hypertension", "prop", "pc4",
  title = "Hypertension: proportionality for PC4", xlab = "Age at onset",
  xlim = c(38,85), xbreaks = seq(40,85,5),
  ylim = c(-0.5,0.5), ybreaks = c(seq(-0.5,0,0.1),seq(0.1,0.5,0.1)),
  text_x = 65, text_y = -0.4, aspect = 4/8)

p_hypertension_prop_pc9 = gof_panel(
  "hypertension", "prop", "pc9",
  title = "Hypertension: proportionality for PC9", xlab = "Age at onset",
  xlim = c(38,85), xbreaks = seq(40,85,5),
  ylim = c(-0.6,0.4), ybreaks = c(seq(-0.6,0,0.1),seq(0.1,0.4,0.1)),
  text_x = 65, text_y = -0.5, aspect = 4/8)

plot = ggarrange(p_hypertension_pc1, p_hypertension_pc4,
                 p_hypertension_pc7, p_hypertension_pc9,
                 p_hypertension_prop_sex, p_hypertension_prop_england,
                 p_hypertension_prop_pc4, p_hypertension_prop_pc9,
                 nrow = 4, ncol = 2,
                 label.y = 0.98, label.x = 0.04,
                 font.label = list(size = 18, color = "black", face = "bold", family = NULL))

ggsave(plot, file = "application/output/FigureS3.pdf", width = 12, height = 14)
ggsave(plot, file = "application/output/FigureS3.eps", width = 12, height = 14)
