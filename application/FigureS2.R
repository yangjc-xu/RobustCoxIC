## ---------------------------------------------------------------------------
## Figure S2: residual processes for the endpoint of diabetes.
##
## Run from the root directory of the reproducibility materials, after
## application/gof_test.R has been run for the two tests shown below:
##   Rscript application/FigureS2.R
##
## Input:  application/output/ff_pc3_diabetes.RData,
##         application/output/prop_sex_diabetes.RData
## Output: application/output/FigureS2.pdf and application/output/FigureS2.eps
## ---------------------------------------------------------------------------

if (!file.exists("application/functions.R")) {
  stop("Please run this script from the root directory of the reproducibility materials.")
}
library(tidyverse)
library(ggpubr)
source("application/functions.R")

p_diabetes_pc3 = gof_panel(
  "diabetes", "ff", "pc3",
  title = "Diabetes: functional form for PC3", xlab = "PC3",
  xlim = c(-8,5), xbreaks = seq(-8,5,1),
  ylim = c(-0.3,0.3), ybreaks = c(seq(-0.3,0,0.1),seq(0.1,0.3,0.1)),
  text_x = 0.2, text_y = -0.24, aspect = 8/8)

p_diabetes_prop_sex = gof_panel(
  "diabetes", "prop", "sex",
  title = "Diabetes: proportionality for sex", xlab = "Age at onset",
  xlim = c(38,85), xbreaks = seq(40,85,5),
  ylim = c(-0.3,0.1), ybreaks = c(seq(-0.3,0,0.05),seq(0.05,0.1,0.05)),
  text_x = 65, text_y = -0.26, aspect = 8/8)

plot = ggarrange(p_diabetes_pc3,p_diabetes_prop_sex,
                 nrow = 1, ncol = 2,
                 label.y = 0.98, label.x = 0.04,
                 font.label = list(size = 18, color = "black", face = "bold", family = NULL))

ggsave(plot, file = "application/output/FigureS2.pdf", width = 10, height = 5)
ggsave(plot, file = "application/output/FigureS2.eps", width = 10, height = 5)
