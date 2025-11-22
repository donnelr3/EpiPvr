rm(list=setdiff(ls(), "filepath"))
setwd("C:/Users/ruair/Desktop/MEE_results_log")

library('EpiPvr')
library(parallel)
library(pbapply)
library(posterior)
#library(dplyr)

nCores <- detectCores() - 1  

CMB=1
CBSI=0


################### FUNCTION ###################################################################################
run_epi_prob_chain <- function(numInsects_vec, al_fits, be_fits, mu_fits,
                               localParams, nChains, virus_run_type) {
  nruns <- nrow(al_fits)
  nCores <- max(1, parallel::detectCores() - 1)
  
  result_vec_byPl <- array(0, dim = c(length(numInsects_vec), nruns, nChains))
  result_vec_byIns <- array(0, dim = c(length(numInsects_vec), nruns, nChains))
  
  cl <- makeCluster(min(nCores, nChains))
  clusterExport(cl, varlist = c("al_fits", "be_fits", "mu_fits", "localParams",
                                "numInsects_vec", "calculate_epidemic_probability", "nruns"),
                envir = environment())
  clusterEvalQ(cl, library(posterior))
  
  # --- progress-enabled parallel call ---
  results_chain_list <- pblapply(1:nChains, function(ggg) {
    result_vec_byPl <- array(0, dim = c(length(numInsects_vec), nruns))
    result_vec_byIns <- array(0, dim = c(length(numInsects_vec), nruns))
    
    pb <- txtProgressBar(min = 0, max = nruns, style = 3)
    cat(sprintf("Starting chain %d...\n", ggg))
    
    logfile <- sprintf("outdata_files/chain_%02d.log", ggg)
    cat(sprintf("Starting chain %d\n", ggg), file = logfile, append = TRUE)
    
    for (ppp in 1:(nruns)) {
      
      if (ppp %% 50 == 0) {
        cat(sprintf("[%s] Chain %d iteration %d/%d\n",
                    format(Sys.time(), "%H:%M:%S"),
                    ggg, ppp, nruns),
            file = logfile, append = TRUE)
      }
      
      setTxtProgressBar(pb, ppp)
      
      
      al_estim <- al_fits[nruns + 1 - ppp, ggg]
      be_estim <- be_fits[nruns + 1 - ppp, ggg]
      mu_estim <- mu_fits[nruns + 1 - ppp, ggg]
      virusParams <- c(al_estim, be_estim, mu_estim)
      
      for (ii in seq_along(numInsects_vec)) {
        numVars <- ((numInsects_vec[ii] + 1) * 3) - 1
        qm_out <- calculate_epidemic_probability(
          numberInsects = numInsects_vec[ii],
          localParameters = localParams,
          virusParameters = virusParams
        )
        result_vec_byPl[ii, ppp] <- qm_out[1]
        result_vec_byIns[ii, ppp] <- qm_out[(numVars - (numInsects_vec[ii] - 1))]
      }
    }
    cat(sprintf("Chain %d complete.\n", ggg), file = logfile, append = TRUE)
    close(pb)
    cat(sprintf("Chain %d completed.\n", ggg))
    
    list(result_vec_byPl = result_vec_byPl,
         result_vec_byIns = result_vec_byIns)
  }, cl = cl)  # <-- progress bar visible here
  stopCluster(cl)
  
  
  for (ggg in 1:nChains) {
    result_vec_byPl[,,ggg]  <- results_chain_list[[ggg]]$result_vec_byPl
    result_vec_byIns[,,ggg] <- results_chain_list[[ggg]]$result_vec_byIns
  }
    
    
    
    
    
  
  # Summarize draws
  data_table_Pl <- array(0, dim = c(3, length(numInsects_vec)))
  data_table_Ins <- array(0, dim = c(3, length(numInsects_vec)))
  target_A1A2 <- array(0, dim = c(2, nruns * nChains, length(numInsects_vec)))
  
  for (ii in seq_along(numInsects_vec)) {
    target_A1 <- matrix(0, nruns, nChains)
    target_A2 <- matrix(0, nruns, nChains)
    for (gg in 1:nChains) {
      target_A1[, gg] <- result_vec_byPl[ii, , gg]
      target_A2[, gg] <- result_vec_byIns[ii, , gg]
    }
    target_A1A2[,,ii] <- rbind(as.vector(target_A1), as.vector(target_A2))
    #ts <- posterior::summarize_draws(posterior::as_draws(t(target_A1A2[,,ii])),quantile_probs = c(0.025, 0.5, 0.975))
    #ts <- posterior::summarize_draws(
    #  posterior::as_draws_matrix(t(target_A1A2[,,ii])),
    #  quantile_probs = c(0.025, 0.5, 0.975)
    #)
    #data_table_Pl[, ii] <- c(ts[1,]$median, ts[1,]$q2.5, ts[1,]$q97.5)
    #data_table_Ins[, ii] <- c(ts[2,]$median, ts[2,]$q2.5, ts[2,]$q97.5)
    draw_mat <- target_A1A2[,,ii]  # 2 x (nruns*nChains)
    data_table_Pl[, ii] <- quantile(draw_mat[1, ], probs = c(0.5, 0.025, 0.975), na.rm = TRUE)
    data_table_Ins[, ii] <- quantile(draw_mat[2, ], probs = c(0.5, 0.025, 0.975), na.rm = TRUE)
  }
  
  # Export results
  dir.create("outdata_files/", showWarnings = FALSE, recursive = TRUE)
  write.table(target_A1A2[1,,], file.path("outdata_files/", paste0("path_", virus_run_type, "_chains_byF_frPlant_ex.dat")))
  write.table(target_A1A2[2,,], file.path("outdata_files/", paste0("path_", virus_run_type, "_chains_byF_frInsect_ex.dat")))
  
  rownames(data_table_Pl) <- rownames(data_table_Ins) <- c("50", "2.5", "97.5")
  colnames(data_table_Pl) <- colnames(data_table_Ins) <- numInsects_vec
  
  write.table(data_table_Pl, file.path("outdata_files/", paste0("path_", virus_run_type, "_summary_byF_frPlant_ex.dat")))
  write.table(data_table_Ins, file.path("outdata_files/", paste0("path_", virus_run_type, "_summary_byF_frInsect_ex.dat")))
  
  invisible(list(
    summary_byPlant = data_table_Pl,
    summary_byInsect = data_table_Ins
  ))
}
################################################################################################################

if (CMB) {                    ######## CMB ANALYSIS ################################
  
  #####
  # 1 #
  #####
  # SIMULATE ASSAY DATA .... following epipvr vignette steps
  ##############################################################################
  set.seed(1000) # CMB
  virusType="PT"
  # set assay structure
  nReps=30 # number of reps
  numWF=10 # number of insects in cohorts
  # virus/insect rates per hr
  alrate=0.1 # acquisition rate
  berate=1 # inoculation rate
  gamrate=0.5 # latency progression rate  (NA if SPT virus)
  murate=0.01 # virus clearance rate
  AAP_lens=c(2,3,3.5,4,4.5,5,6,7,8); # variable duration vectors, hours AAP sub-assay
  LAP_lens=c(0.5,1,2,3,4,5,6,7,8);                                        # LAP sub-assay
  IAP_lens=c(5,10,15,20,25,30,40,50,60)/60;                               # IAP sub-assay
  
  AAP_Reps=rep(nReps,length(AAP_lens));  # vectors for number of replicates
  LAP_Reps=rep(nReps,length(LAP_lens));
  IAP_Reps=rep(nReps,length(IAP_lens));
  
  AAP_Infs_zeros=rep(0,length(AAP_lens));  # vectors for num infected test plants-initially 0
  LAP_Infs_zeros=rep(0,length(LAP_lens));
  IAP_Infs_zeros=rep(0,length(IAP_lens));
  
  AAPinput=rbind(AAP_lens, AAP_Reps, AAP_Infs_zeros)  # group structural vectors for input
  LAPinput=rbind(LAP_lens, LAP_Reps, LAP_Infs_zeros)
  IAPinput=rbind(IAP_lens, IAP_Reps, IAP_Infs_zeros)
  
  # default durations of acquisition, latent and inoculation periods
  T_A=4
  T_L=2
  T_I=6
  
  # Place '-1' in the varied assay component
  AAPfixedComponent=c(-1,T_L,T_I)
  LAPfixedComponent=c(T_A,-1,T_I)
  IAPfixedComponent=c(T_A,T_L,-1)
  ddur_mat=rbind(AAPfixedComponent,LAPfixedComponent,IAPfixedComponent)
  
  # simulate and hence populate the number of infected test plants in the AAP subassay
  assay1=AP_assay_simulator(AAPinput,
                            AAPfixedComponent,
                            numWF, c(alrate,berate,gamrate,murate),isVerbose=0,'PT')  
  # simulate and hence populate the number of infected test plants in the LAP subassay
  assay2=AP_assay_simulator(LAPinput,
                            LAPfixedComponent,
                            numWF, c(alrate,berate,gamrate,murate),isVerbose=0,'PT') 
  # simulate and hence populate the number of infected test plants in the IAP subassay
  assay3=AP_assay_simulator(IAPinput,
                            IAPfixedComponent,
                            numWF, c(alrate,berate,gamrate,murate),isVerbose=0,'PT') 
  
  ap_data_sim_PT=list(d_AAP=assay1,      # AAP sub-assay structure 
                   d_LAP=assay2,      # LAP sub-assay structure  (OMIT IF SPT VIRUS) 
                   d_IAP=assay3,      # IAP sub-assay structure 
                   d_durations=ddur_mat,  # fixed durations
                   d_vectorspp=numWF,     # insects in a cohort
                   d_virusType=virusType) # virus code - element of {'SPT','PT'}
  
  # Adding virus parameters as attributes nb. -1 indicates unknown virus parameters
  attr(ap_data_sim_PT, "alpha") <- alrate # acquisition rate
  attr(ap_data_sim_PT, "beta") <- berate  # inocualation rate
  if (!is.na(gamrate)) {
    attr(ap_data_sim_PT, "gamma") <- gamrate  # latent progression rate(OMIT IF SPT VIRUS) 
  }
  attr(ap_data_sim_PT, "mu") <- murate  # insect recovery rate rate
  
  print(ap_data_sim_PT)
  #save to rda file
  save(ap_data_sim_PT, file = "outdata_files/ap_data_sim_PT.rda")
  ##############################################################################
  
  #####
  # 2 #
  #####
  # ANALYSE ASSAY DATA .... 
  ##############################################################################
  set.seed(1001) # CMB
  
  nChains=4
  lsEst_in=40;
  numPtsPdin=1;
  numParams=4;
  
  # first run on simulated dataset - epi parameter analysis only
  mcmcOptions=c(1000,2000)
  EVPT_sim=estimate_virus_parameters_PT(ap_data_sim_PT,lsEst_in,numPtsPdin,mcmcOptions,numChainsIn = 4,mc.parallel=1)
  saveRDS(EVPT_sim, file="outdata_files/EVPT_sim.rds") # save the chains
  
  assign("last.warning", NULL, envir = baseenv())
  
  # then run on the published cmb dataset - epi parameter analysis AND epidemic risk calculation
  mcmcOptions=c(1000,2000)
  ap_data_PTpub <- get(load("data/PT_PUB_AccessExp.rda"))
  ap_data_PTpub$d_virusType="PT"
  ap_data_PTpub$d_LAP=0*ap_data_PTpub$d_LAP   # WE DO THIS BECAUSE DUBERN 94 used non-cassava host plant for intermediate sub-assay
  EVPT_pub=estimate_virus_parameters_PT(ap_data_PTpub,lsEst_in,numPtsPdin,mcmcOptions,numChainsIn = 4,mc.parallel=1)
  saveRDS(EVPT_pub, file="outdata_files/EVPT_pub.rds") # save the chains
  target=EVPT_pub$array1
  
  assign("last.warning", NULL, envir = baseenv())
 
  # epidemic risk calculation 
  # P EPIDEMIC inferences  # accessing the chains from estimate_virus_parameters  # VIRUS PARAMETERS
  al_fits=target[, , "al[1]",drop=TRUE]*24 # convert from per min to per day
  be_fits=target[, , "be[1]",drop=TRUE]*24
  mu_fits=target[, , "mu[1]",drop=TRUE]*24
  print(sum((al_fits<0)))
  print(sum((be_fits<0)))
  print(sum((mu_fits<0)))
  # LOCAL PARAMETERS  # set the local parameters (per day)
  thet_external <- 0.45 # dispersal
  r_external  <- 1/28 # roguing
  bf_external <- 1/14  # vector mortality
  h_external <- 1/365  # harvesting rate
  nu_pl_external <- 1/14  # plant latent progression rate (1/plant latent period)
  localParams=c(thet_external, r_external, h_external, bf_external, nu_pl_external)
  # generate p Epdemic results for various insect burdens (per plant)
  numInsects_vec_cmb=c(1,2,3)
  # epidemic probability calculator on the entire chain
  results=run_epi_prob_chain(numInsects_vec_cmb, al_fits, be_fits, mu_fits, localParams, nChains, "cmb")

}else if (CBSI) {                    ######## CBSI ANALYSIS ################################
  
  #####
  # 1 #
  #####
  # SIMULATE ASSAY DATA .... following epipvr vignette steps
  ##############################################################################
  set.seed(1011) # CBSI
  virusType="SPT"
  # set assay structure
  nReps=30 # number of reps
  numWF=20 # number of insects in cohorts
  # virus/insect rates per hr
  alrate=0.1 # acquisition rate
  berate=1 # inoculation rate
  gamrate=NA # latency progression rate  (NA if SPT virus)
  murate=1 # virus clearance rate
  AAP_lens=c(2,3,3.5,4,4.5,5,6,7,8); # variable duration vectors, hours AAP sub-assay
  IAP_lens=c(5,10,15,20,25,30,40,50,60)/60;                               # IAP sub-assay
  
  AAP_Reps=rep(nReps,length(AAP_lens));  # vectors for number of replicates
  IAP_Reps=rep(nReps,length(IAP_lens));
  
  AAP_Infs_zeros=rep(0,length(AAP_lens));  # vectors for num infected test plants-initially 0
  IAP_Infs_zeros=rep(0,length(IAP_lens));
  
  AAPinput=rbind(AAP_lens, AAP_Reps, AAP_Infs_zeros)  # group structural vectors for input
  IAPinput=rbind(IAP_lens, IAP_Reps, IAP_Infs_zeros)
  
  # default durations of acquisition, latent and inoculation periods
  T_A=4
  T_I=6
  
  # Place '-1' in the varied assay component
  AAPfixedComponent=c(-1,0,T_I)
  IAPfixedComponent=c(T_A,0,-1)
  
  # simulate and hence populate the number of infected test plants in the AAP subassay
  assay1=AP_assay_simulator(AAPinput,
                            AAPfixedComponent,
                            numWF, c(alrate,berate,gamrate,murate),isVerbose=0,'SPT')  
 
  # simulate and hence populate the number of infected test plants in the IAP subassay
  assay3=AP_assay_simulator(IAPinput,
                            IAPfixedComponent,
                            numWF, c(alrate,berate,gamrate,murate),isVerbose=0,'SPT') 
  
  ddur_mat=rbind(AAPfixedComponent[c(1,3)],IAPfixedComponent[c(1,3)])
  
  ap_data_sim_SPT=list(d_AAP=assay1,      # AAP sub-assay structure 
                   d_IAP=assay3,      # IAP sub-assay structure 
                   d_durations=ddur_mat,  # fixed durations
                   d_vectorspp=numWF,     # insects in a cohort
                   d_virusType=virusType) # virus code - element of {'SPT','PT'}
  
  # Adding virus parameters as attributes nb. -1 indicates unknown virus parameters
  # note that dataset format does not include gamma atttibute for SPT virus
  attr(ap_data_sim_SPT, "alpha") <- alrate 
  attr(ap_data_sim_SPT, "beta") <- berate 
  attr(ap_data_sim_SPT, "mu") <- murate 
  
  print(ap_data_sim_SPT)
  #save to rda file
  save(ap_data_sim_SPT, file = "outdata_files/ap_data_sim_SPT.rda")
  ##############################################################################
  
  
  
  #####
  # 2 #
  #####
  # ANALYSE ASSAY DATA .... 
  ##############################################################################
  set.seed(1002) # CBSI
  
  nChains=4
  lsEst_in=40;
  numPtsPdin=1;
  numParams=3;
  
  # first run on simulated dataset - epi parameter analysis only
  mcmcOptions=c(4500,6000)
  EVSPT_sim=estimate_virus_parameters_SPT(ap_data_sim_SPT,lsEst_in,numPtsPdin,mcmcOptions,numChainsIn = 4,mc.parallel=1)
  saveRDS(EVSPT_sim, file="outdata_files/EVSPT_sim") # save the chains
  
  # then run on the published cmb dataset - epi parameter analysis AND epidemic risk calculation
  mcmcOptions=c(4500,6000)
  ap_data_SPTpub <- get(load("data/SPT_PUB_AccessExp.rda"))
  ap_data_SPTpub$d_virusType="SPT"
  EVSPT_pub=estimate_virus_parameters_SPT(ap_data_SPTpub,lsEst_in,numPtsPdin,mcmcOptions,numChainsIn = 4,mc.parallel=1)
  saveRDS(EVSPT_pub, file="outdata_files/EVSPT_pub") # save the chains
  target=EVSPT_pub$array1
  

  
  #### INFERENCES #######
  # P EPIDEMIC inferences # accessing the chains from estimate_virus_parameters # VIRUS PARAMETERS
  al_fits=target[, , "al[1]",drop=TRUE]*24 # convert from per min to per day 
  be_fits=target[, , "be[1]",drop=TRUE]*24
  mu_fits=target[, , "mu[1]",drop=TRUE]*24
  print(sum((al_fits<0)))
  print(sum((be_fits<0)))
  print(sum((mu_fits<0)))
  # LOCAL PARAMETERS
  # set the local parameters (per day)
  thet_external <- 0.45 # dispersal
  r_external  <- 1/28 # roguing
  bf_external <- 1/14  # vector mortality
  h_external <- 1/365  # harvesting rate
  nu_pl_external <- 1/14  # plant latent progression rate (1/plant latent period)
  localParams=c(thet_external, r_external, h_external, bf_external, nu_pl_external)
  # generate p Epdemic results for various insect burdens (per plant)
  numInsects_vec_cbsi=c(4,7,10)
  # epidemic probability calculator on the entire chain
  results=run_epi_prob_chain(numInsects_vec_cbsi, al_fits, be_fits, mu_fits, localParams, nChains, "cbsi")
}

