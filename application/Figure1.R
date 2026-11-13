## ---------------------------------------------------------------------------
## Figure 1: standard error estimates and p-values of the model-based and
## robust score statistics for diabetes and hypertension. The script also
## prints the numbers of significant genetic variants reported in Section 7.
##
## Run from the root directory of the reproducibility materials, after
## application/score_test.R has been run for both endpoints:
##   Rscript application/Figure1.R
##
## Input:  application/output/score_test_{endpoint}.RData
## Output: application/output/Figure1.pdf and application/output/Figure1.eps
## ---------------------------------------------------------------------------

if (!file.exists("application/functions.R")) {
  stop("Please run this script from the root directory of the reproducibility materials.")
}
library(tidyverse)
library(latex2exp)
library(ggpubr)

alpha1 = 0.05/220
alpha2 = 5e-8

## Two panels for one endpoint: standard error estimates and p-values.
## Column 2 of variance and p_value is model-based, and column 3 is robust.
score_panels = function(endpoint, title, se_limits, se_breaks){
  # loads statistic, variance, p_value, and variant_name
  load(sprintf("application/output/score_test_%s.RData", endpoint))
  
  for (level in c(alpha1, alpha2)) {
    cat(title, ", significance level", level, ": model-based score test,",
        length(which(p_value[,2]<level)), "significant variants; robust score test,",
        length(which(p_value[,3]<level)), "significant variants\n")
  }
  
  df_se = data.frame(variant_name = variant_name, se_naive = sqrt(variance[,2]), se_robust = sqrt(variance[,3]))
  
  p_se = ggplot(data = df_se) + 
    geom_point(aes(x = se_robust, y = se_naive)) +
    labs(title = title) +
    scale_x_continuous(limits = se_limits, breaks = se_breaks, name = "Robust standard error estimates") +
    scale_y_continuous(limits = se_limits, breaks = se_breaks, name = "Model-based standard error estimates")+
    geom_abline(slope = 1, intercept = 0, colour = "#E41A1C", linewidth = 0.8) +
    theme_light() +
    theme(legend.key = element_rect(colour = NA, fill = NA),
          legend.title = element_blank(),
          legend.background = element_blank(), 
          legend.key.height = unit(10, "pt"),
          legend.key.width = unit(20, "pt"),
          axis.title.x = element_text(size = 15, margin = margin(t=8)),
          axis.text = element_text(size = 15),
          legend.text = element_text(size = 10),
          plot.margin = margin(t=1,1,1,1, "lines"),
          aspect.ratio=10/10,
          title = element_text(size = 15))
  
  df_pval = data.frame(variant_name = variant_name, p_naive = p_value[,2], p_robust = p_value[,3]) %>%
    filter(p_naive>0, p_robust>0)
  
  p_val = ggplot(data = df_pval) + 
    geom_point(aes(x = -log10(p_robust), y = -log10(p_naive))) +
    labs(title = title) +
    scale_x_continuous(limits = c(0,16), breaks =  seq(0,16,2), name = TeX(r'(Robust $-log_{10}$($p$-values))')) +
    scale_y_continuous(limits = c(0,16), breaks = seq(0,16,2), name = TeX(r'(Model-based $-log_{10}$($p$-values))'))+
    geom_abline(slope = 1, intercept = 0, colour = "#E41A1C", linewidth = 0.8) +
    geom_hline(yintercept = -log10(alpha2), color = "#377EB8", linetype="dashed", linewidth = 0.8)+
    geom_vline(xintercept = -log10(alpha2), color = "#377EB8", linetype="dashed", linewidth = 0.8)+
    geom_hline(yintercept = -log10(alpha1), color = "#4DAF4A", linetype="dashed", linewidth = 0.8)+
    geom_vline(xintercept = -log10(alpha1), color = "#4DAF4A", linetype="dashed", linewidth = 0.8)+
    theme_light() +
    theme(legend.key = element_rect(colour = NA, fill = NA),
          legend.title = element_blank(),
          legend.background = element_blank(), 
          legend.key.height = unit(10, "pt"),
          legend.key.width = unit(20, "pt"),
          axis.title.x = element_text(size = 15, margin = margin(t=8)),
          axis.text = element_text(size = 15),
          legend.text = element_text(size = 10),
          plot.margin = margin(t=1,1,1,1, "lines"),
          aspect.ratio=10/10,
          title = element_text(size = 15))
  
  list(p_se = p_se, p_val = p_val)
}

panels_dibt = score_panels("diabetes", title = "Diabetes",
                           se_limits = c(0.03,0.18), se_breaks = seq(0.03,0.18,0.03))
panels_hypt = score_panels("hypertension", title = "Hypertension",
                           se_limits = c(0.25,0.55), se_breaks = seq(0.25,0.55,0.05))

p = ggarrange(panels_dibt$p_se, panels_dibt$p_val, 
              panels_hypt$p_se, panels_hypt$p_val, 
              ncol = 2, nrow = 2)
ggsave(p, file = "application/output/Figure1.pdf", width = 10, height = 10, units = "in")
ggsave(p, file = "application/output/Figure1.eps", width = 10, height = 10, units = "in")
