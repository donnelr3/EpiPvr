rm(list=setdiff(ls(), "filepath"))
setwd("C:/Users/ruair/Desktop/MEE_results_log")

#library('EpiPvr')
#library(parallel)
#library(pbapply)
library(posterior)
library(dplyr)


#--------------------------------------------------------------------------#
#--------------------------------------------------------------------------#
#--------------------------------------------------------------------------#


# PT

EVPT_pub_in=readRDS(file="outdata_files/EVPT_pub.rds") # load the chains
draws_df_PT_pub <- as_draws_df(EVPT_pub_in$array1)
print(draws_df_PT_pub)
summary_stats_PT_pub <- summarize_draws(
  draws_df_PT_pub,
  ~quantile2(.x, probs = c(0.025, 0.5, 0.975)),
  ~mean(.x),
  ~sd(.x)
)
ci_parameter_table_PT_pub = summary_stats_PT_pub[,-6]

EVPT_sim_in=readRDS(file="outdata_files/EVPT_sim.rds") # load the chains
draws_df_PT_sim <- as_draws_df(EVPT_sim_in$array1)
print(draws_df_PT_sim)
summary_stats_PT_sim <- summarize_draws(
  draws_df_PT_sim,
  ~quantile2(.x, probs = c(0.025, 0.5, 0.975)),
  ~mean(.x),
  ~sd(.x)
)
ci_parameter_table_PT_sim = summary_stats_PT_sim[,-6]

# SPT

EVSPT_pub_in=readRDS(file="outdata_files/EVSPT_pub") # load the chains
draws_df_SPT_pub <- as_draws_df(EVSPT_pub_in$array1)
print(draws_df_SPT_pub)
summary_stats_SPT_pub <- summarize_draws(
  draws_df_SPT_pub,
  ~quantile2(.x, probs = c(0.025, 0.5, 0.975)),
  ~mean(.x),
  ~sd(.x)
)
ci_parameter_table_SPT_pub = summary_stats_SPT_pub[,-6]

EVSPT_sim_in=readRDS(file="outdata_files/EVSPT_sim") # load the chains
draws_df_SPT_sim <- as_draws_df(EVSPT_sim_in$array1)
print(draws_df_SPT_sim)
summary_stats_SPT_sim <- summarize_draws(
  draws_df_SPT_sim,
  ~quantile2(.x, probs = c(0.025, 0.5, 0.975)),
  ~mean(.x),
  ~sd(.x)
)
ci_parameter_table_SPT_sim = summary_stats_SPT_sim[,-6]

################################################################################

ins_EVPT_pub_in=read.table("outdata_files/path_cmb_chains_byF_frInsect.dat")
pl_EVPT_pub_in=read.table("outdata_files/path_cmb_chains_byF_frPlant.dat")
n_chain <- 4
n_iter <- nrow(ins_EVPT_pub_in)/n_chain
n_param <- ncol(ins_EVPT_pub_in)
ins_draws_array_PT <- array(
  data = as.matrix(ins_EVPT_pub_in),
  dim = c(n_iter, n_chain, n_param),
  dimnames = list(
    iteration = NULL,
    chain = paste0("chain", 1:n_chain),
    variable = colnames(ins_EVPT_pub_in)
  )
)
pl_draws_array_PT <- array(
  data = as.matrix(pl_EVPT_pub_in),
  dim = c(n_iter, n_chain, n_param),
  dimnames = list(
    iteration = NULL,
    chain = paste0("chain", 1:n_chain),
    variable = colnames(pl_EVPT_pub_in)
  )
)
ins_draws_df_PT <- as_draws_array(ins_draws_array_PT)
pl_draws_df_PT <- as_draws_array(pl_draws_array_PT)
ins_summary_stats_PT <- summarize_draws(
  ins_draws_df_PT,
  ~quantile2(.x, probs = c(0.025, 0.5, 0.975)),
  ~mean(.x),
  ~sd(.x)
)
pl_summary_stats_PT <- summarize_draws(
  pl_draws_df_PT,
  ~quantile2(.x, probs = c(0.025, 0.5, 0.975)),
  ~mean(.x),
  ~sd(.x)
)
ins_ci_parameter_table_PT = ins_summary_stats_PT[,-6]
pl_ci_parameter_table_PT = pl_summary_stats_PT[,-6]



ins_EVSPT_pub_in=read.table("outdata_files/path_cbsi_chains_byF_frInsect.dat")
pl_EVSPT_pub_in=read.table("outdata_files/path_cbsi_chains_byF_frPlant.dat")
n_chain <- 4
n_iter <- nrow(ins_EVSPT_pub_in)/n_chain
n_param <- ncol(ins_EVSPT_pub_in)
ins_draws_array_SPT <- array(
  data = as.matrix(ins_EVSPT_pub_in),
  dim = c(n_iter, n_chain, n_param),
  dimnames = list(
    iteration = NULL,
    chain = paste0("chain", 1:n_chain),
    variable = colnames(ins_EVPT_pub_in)
  )
)
pl_draws_array_SPT <- array(
  data = as.matrix(pl_EVSPT_pub_in),
  dim = c(n_iter, n_chain, n_param),
  dimnames = list(
    iteration = NULL,
    chain = paste0("chain", 1:n_chain),
    variable = colnames(pl_EVSPT_pub_in)
  )
)
ins_draws_df_SPT <- as_draws_array(ins_draws_array_SPT)
pl_draws_df_SPT <- as_draws_array(pl_draws_array_SPT)
ins_summary_stats_SPT <- summarize_draws(
  ins_draws_df_SPT,
  ~quantile2(.x, probs = c(0.025, 0.5, 0.975)),
  ~mean(.x),
  ~sd(.x)
)
pl_summary_stats_SPT <- summarize_draws(
  pl_draws_df_SPT,
  ~quantile2(.x, probs = c(0.025, 0.5, 0.975)),
  ~mean(.x),
  ~sd(.x)
)
ins_ci_parameter_table_SPT = ins_summary_stats_SPT[,-6]
pl_ci_parameter_table_SPT = pl_summary_stats_SPT[,-6]











#### validation plots #### validation plots #### validation plots #### validation plots #### validation plots

# =========================================================
#  SPT SIMULATED  --- pair of plots (AAP and IAP)
# =========================================================

pdf("APdata_valdn_plot_SPT_SIMULATED.pdf")
par(mfrow = c(1, 2), mar = c(5, 4, 4, 0))  # first plot wider left margin

# ---- AAP ----
plot(EVSPT_sim_in$array3$lenA, EVSPT_sim_in$array3$propA,
     xlab = "AAP duration, hours",
     ylab = "Prop. test plants SPT infected",
     ylim = c(0, 1),
     pch = 1)
lines(EVSPT_sim_in$array3$lenA, EVSPT_sim_in$array3$simulLA)
lines(EVSPT_sim_in$array3$lenA, EVSPT_sim_in$array3$simulUA)
legend(EVSPT_sim_in$array3$lenA[1], 0.5,
       legend = c("Experiment", "95% Cr.I."),
       pch = c(1, NA),
       lty = c(NA, 1),
       col = "black",
       cex = 0.8,
       bty = "n")

# ---- IAP ----
par(mar = c(5, 2, 4, 2))  # smaller left margin
plot(EVSPT_sim_in$array4$lenI, EVSPT_sim_in$array4$propI,
     xlab = "IAP duration, hours",
     ylab = "Prop. test plants SPT infected",
     ylim = c(0, 1),
     pch = 1)
lines(EVSPT_sim_in$array4$lenI, EVSPT_sim_in$array4$simulLI)
lines(EVSPT_sim_in$array4$lenI, EVSPT_sim_in$array4$simulUI)
legend(EVSPT_sim_in$array4$lenI[1], 0.5,
       legend = c("Experiment", "95% Cr.I."),
       pch = c(1, NA),
       lty = c(NA, 1),
       col = "black",
       cex = 0.8,
       bty = "n")

dev.off()


# =========================================================
#  PT SIMULATED  --- triple plots (AAP, LAP, IAP)
#      2 columns top + 1 centered below (equal width)
# =========================================================

pdf("APdata_valdn_plot_PT_SIMULATED.pdf")

# layout matrix with a blank cell in bottom-left, third plot in bottom-right
layout(matrix(c(1, 2,
                0, 3), ncol = 2, byrow = TRUE),
       heights = c(1, 1), widths = c(1, 1))

# --- AAP (top-left) ---
par(mar = c(5, 4, 4, 0))
plot(EVPT_sim_in$array3$lenA, EVPT_sim_in$array3$propA,
     xlab = "AAP duration, hours",
     ylab = "Prop. test plants PT infected",
     ylim = c(0, 1),
     pch = 1)
lines(EVPT_sim_in$array3$lenA, EVPT_sim_in$array3$simulLA)
lines(EVPT_sim_in$array3$lenA, EVPT_sim_in$array3$simulUA)

# --- LAP (top-right) ---
par(mar = c(5, 4, 4, 0))
plot(EVPT_sim_in$array4$lenL, EVPT_sim_in$array4$propL,
     xlab = "LAP duration, hours",
     ylab = "Prop. test plants PT infected",
     ylim = c(0, 1),
     pch = 1)
lines(EVPT_sim_in$array4$lenL, EVPT_sim_in$array4$simulLL)
lines(EVPT_sim_in$array4$lenL, EVPT_sim_in$array4$simulUL)

# --- IAP (bottom-right only, same width as above) ---
par(mar = c(5, 4, 4, 0))
plot(EVPT_sim_in$array5$lenI, EVPT_sim_in$array5$propI,
     xlab = "IAP duration, hours",
     ylab = "Prop. test plants PT infected",
     ylim = c(0, 1),
     pch = 1)
lines(EVPT_sim_in$array5$lenI, EVPT_sim_in$array5$simulLI)
lines(EVPT_sim_in$array5$lenI, EVPT_sim_in$array5$simulUI)

legend("bottomright",
       legend = c("Experiment", "95% Cr.I."),
       pch = c(1, NA),
       lty = c(NA, 1),
       col = "black",
       cex = 0.8,
       bty = "n")

dev.off()
