## ---------------------------------------------------------------------------
## Shared settings and functions for the UK Biobank data application
## (Figure 1, Figures S2 and S3). This file only defines objects and functions;
## it is sourced by the other scripts in application/.
## ---------------------------------------------------------------------------

## Endpoint-specific settings.
##   data_file:     location of the analysis data file (see DATA_DICTIONARY.md).
##   variant_cols:  columns of the data file holding the 220 genetic variants.
##   round_divisor: divisor applied to the interval endpoints before rounding
##                  them down to whole years (used for the proportional
##                  hazards checks); diabetes.tab records ages in days and
##                  hypertension.tab records ages in years.
endpoint_config = list(
  diabetes = list(data_file = "data/diabetes.tab",
                  variant_cols = 17:236,
                  round_divisor = 365),
  hypertension = list(data_file = "data/hypertension.tab",
                      variant_cols = 4:223,
                      round_divisor = 1)
)

## Baseline covariates, in the order of the columns of the design matrix X of
## the null model.
covariate_names = c("sex", "white", "england", paste0("pc", 1:10))

## Covariates examined in the goodness-of-fit tests of each null model.
##   ff:   functional-form checks, with the range of the grid of 101 cut-off
##         values for each covariate.
##   prop: proportional hazards checks.
gof_config = list(
  diabetes = list(
    ff = list(pc1 = c(-20,-5), pc3 = c(-8,5), pc7 = c(-10,10)),
    prop = c("sex", "pc1", "pc3", "pc7")
  ),
  hypertension = list(
    ff = list(pc1 = c(-20,-5), pc4 = c(-10,15), pc6 = c(-5,5),
              pc7 = c(-10,10), pc8 = c(-10,10), pc9 = c(-15,10)),
    prop = c("sex", "white", "england", "pc1", "pc4", "pc6", "pc7", "pc8", "pc9")
  )
)

check_endpoint = function(endpoint){
  if (is.na(endpoint) || !(endpoint %in% names(endpoint_config))) {
    stop("endpoint must be 'diabetes' or 'hypertension'")
  }
  endpoint
}

## One panel of Figures S2 and S3: the observed residual process and 20
## simulated processes saved by gof_test.R.
##   sim_rows: rows of W_simu1 (simulated processes) to display.
gof_panel = function(endpoint, type, covariate, title, xlab,
                     xlim, xbreaks, ylim, ybreaks, text_x, text_y,
                     aspect, sim_rows = 1:20){
  # loads W, W_simu1, x, and pval
  load(sprintf("application/output/%s_%s_%s.RData", type, covariate, endpoint))
  
  df_observed = data.frame(x = x, W = W)
  df_simulated = data.frame(
    x = rep(x, times = length(sim_rows)),
    W = as.vector(t(W_simu1[sim_rows, , drop = FALSE])),
    curve = rep(seq_along(sim_rows), each = length(x))
  )
  pval_label = if (pval < 0.001) "P-value < 0.001" else sprintf("P-value = %.3f", pval)
  
  ggplot() +
    theme_light() +
    geom_line(data = df_simulated,
              aes(x = x, y = W, group = curve, color = "Simulated", linewidth = "Simulated", linetype = "Simulated")) +
    geom_line(data = df_observed,
              aes(x = x, y = W, color = "Observed", linewidth = "Observed", linetype = "Observed")) +
    labs(x = xlab, y = "Cumulative statistic", title = title) +
    scale_x_continuous(limits = xlim, breaks = xbreaks, expand = expansion(0,0)) +
    scale_y_continuous(limits = ylim, breaks = ybreaks, expand = expansion(0,0)) +
    scale_color_manual("", values = c(alpha(c("orange", "blue"), 1)), breaks = c("Simulated", "Observed")) +
    scale_linetype_manual("", values = c("dashed", "solid"), breaks = c("Simulated", "Observed")) +
    scale_linewidth_manual("", values = c(0.5,1), breaks = c("Simulated", "Observed")) +
    annotate("text", x = text_x, y = text_y, label = pval_label, hjust = 0, size = 4.2, color = "red") +
    theme(legend.position = "none",
          legend.key = element_rect(colour = NA, fill = NA),
          legend.background = element_blank(), 
          legend.key.height = unit(17, "pt"),
          legend.key.width = unit(30, "pt"),
          axis.title.x = element_text(margin = margin(t=8)),
          axis.text.x = element_text(angle = 0, vjust = 0.6),
          axis.text = element_text(size = 12),
          plot.margin = margin(t=0.5,0.5,0.5,0.5, "lines"),
          aspect.ratio = aspect,
          title = element_text(size = 15))
}
