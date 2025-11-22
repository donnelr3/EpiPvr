#' Estimate PT virus parameters using Bayesian inference with rstan
#'
#' Runs MCMC sampling using a precompiled Stan model for PT plant virus. Analyses PT access period data estimating 5 parameters: (\code{alpha[1]}, \code{beta[1]},
#' \code{gamma[1]},\code{mu[1]},\code{bd[1]} - acquisition, inoculation, latent progression, insect recovery and insect lab mortality rates, respectively).
#'
#' @param data A list containing the AP experiment data for Stan (required).
#' @param D_LSin Upper bound on insect vector lifespan in days sets the vector survival discretization (optional, default = 50).
#' @param D_numPtsPdin Number of points per day of insect vector sets the vector survival discretization (optional, default = 1).
#' @param mcmcOptions A numeric vector of length 2:
#'   The first element specifies the number of warm-up iterations (optional, default = 500),
#'   and the second element specifies the total number of iterations (optional, default = 1000).
#' @param numChainsIn Numeric: number of Markov chains (optional, default = 4).
#' @param mc.parallel Binary: whether to use parallelisation for sampling (optional, default = 0, i.e. 1 core only).
#' @return A list with seven elements:
#'   \describe{
#'     \item{array0}{MCMC chains for estimated parameters (\code{alpha[1]}, \code{beta[1]}, \code{gamma[1]/mu[1]},\code{bd[1]} - and additionally \code{lp__} for reference only).}
#'     \item{array1}{MCMC chains for estimated parameters (\code{alpha[1]}, \code{beta[1]}, \code{gamma[1]},\code{mu[1]}).}
#'     \item{array2}{summary_stats (rhat, ess_bulk, ess_tail). \code{[1]}: R-hat statistic for convergence (should be close to 1); \code{[2-3]}: statistics for n-eff.}
#'     \item{array3}{Validation dataset: AAP sub-assay input test plant data, forward simulated 2.5th percentile, forward simulated 97.5th percentile.}
#'     \item{array4}{Validation dataset: LAP sub-assay input test plant data, forward simulated 2.5th percentile, forward simulated 97.5th percentile.}
#'     \item{array5}{Validation dataset: IAP sub-assay input test plant data, forward simulated 2.5th percentile, forward simulated 97.5th percentile.}
#'     \item{array6}{Mean and sd Bayesian R-squared values for model fit assessment, for AAP, LAP and IAP assays.}
#'     \item{converge_results}{A list containing 2 sampler diagnostics, the number of overall divergent transitions and if maximum tree depth has been exceeded.
#'                      See also screen print for acceptability of energy Bayesian fraction of missing information (E-BFMI).}
#'   }
#'
#' @details
#' **Interpreting model output - check stability first**
#'
#' Warnings will be printed to the screen if there are issues with model fitting.
#' **Do not suppress warnings** (e.g., via \code{suppressWarnings()}), as they contain essential information about
#' potential convergence problems.
#'
#' Before interpreting any model results, always check the following **core diagnostics**:
#'
#' 1. **R-hat** - Should be close to 1 for all parameters; larger values indicate non-convergence
#'    (\href{https://mc-stan.org/learn-stan/diagnostics-warnings.html}{Stan R-hat documentation}).
#' 2. **Effective Sample Size (n_eff)** - Very low values suggest high autocorrelation or insufficient sampling
#'    (\href{https://mc-stan.org/docs/2_19/reference-manual/effective-sample-size-section.html}{Stan ESS documentation}).
#' 3. **Divergent transitions** - Count should be 0; any non-zero count requires investigation
#'    (\href{https://mc-stan.org/docs/2_19/reference-manual/divergent-transitions}{Stan divergences documentation}).
#' 4. **Forward Simulation from Fitted Model (simulated_data)** - Presented conveniently for plotting to assess posterior predictive fit.
#'    - \code{array3}: AAP forward simulation 95% credible intervals and original data
#'    - \code{array4}: LAP forward simulation 95% credible intervals and original data
#'    - \code{array5}: IAP forward simulation 95% credible intervals and original data
#'    (\href{https://mc-stan.org/docs/stan-users-guide/posterior-predictive-checks.html}{Stan posterior predictive checks}).
#'
#' If any of these core diagnostics fail, the model fit may not be trustworthy.
#'
#'
#' **Additional diagnostics for troubleshooting**
#'
#' These are useful when core checks indicate problems, or for deeper exploration of model behaviour:
#'
#' - **Bayesian R^2** (\code{r2_bayesA}, \code{r2_bayesL}, \code{r2_bayesI}) - Measures explanatory power;
#'   low values suggest poor fit to the data (i.e., 3 values for AAP assay, LAP assay,IAP assay).
#' - **Max treedepth exceeded** - Number of times the sampler hit the maximum tree depth; should be 0
#'   (\href{https://mc-stan.org/learn-stan/diagnostics-warnings.html}{Stan max treedepth documentation}).
#' - **EBFMI** (Energy Bayesian Fraction of Missing Information) - Low values indicate poor exploration
#'   of the posterior.
#'
#' See the package vignette for worked examples of checking and interpreting these diagnostics.
#'
#' **Preprint reference**
#'
#' The models implemented in this function follow Donnelly et al. (2025, preprint), originally implemented in the EpiPv GitHub package.
#'
#' @importFrom rstan sampling summary get_sampler_params check_energy
#' @importFrom posterior as_draws_array subset_draws summarise_draws
#' @importFrom stats sd var
#' @references
#' Donnelly, R., Tankam, I. & Gilligan, C. (2025).
#' "Plant pathogen profiling with the EpiPv package."
#' EcoEvoRxiv, \doi{10.32942/X29K9P}.
#'
#' When available, please cite the published version.
#'
#' @examples
#' data("ap_data_sim_PT", package = "EpiPvr")
#' # run low warm-up and iterations (mcmcOptions) for quick example only
#' EVPT_pub=estimate_virus_parameters_PT(
#'   data=ap_data_sim_PT,D_LSin=5,D_numPtsPdin=1,
#'   mcmcOptions=c(25,50),numChainsIn=1,mc.parallel=0
#' )
#' @export
estimate_virus_parameters_PT <- function(data,D_LSin=50,D_numPtsPdin=1,mcmcOptions=c(500,1000),numChainsIn=4,mc.parallel=0){



  numWarm=mcmcOptions[1]
  numIter=mcmcOptions[2]
  adaptVal=0.95 #1-(10^-3)
  treeDepth=10 #15

  numChains = numChainsIn

  message(paste('numChains =',numChainsIn))

  ##### diags should be negative in assay durations matrix (corresponds to separate duration vector that is varied) #####
  if (sum(diag(data$d_durations)>0)!=0){
    stop('input error! see data d_durations');
  }

  ##### data should only involve SPT virus #####
  if (!identical(data$d_virusType,"PT")){
    stop('input error! see data d_virusType');
  }

  ##### data must have acquisition access and inoculation access to be valid AP data #####
  if (any((data$d_durations[1,3]==0),(data$d_durations[3,1]==0))){
    stop('input error! see data d_durations');
  }



  # assembling all the data inputs for the stan estimation
  dat3= list(D_NumGrps=c(dim(data$d_AAP)[2],dim(data$d_LAP)[2],dim(data$d_IAP)[2]),#length(AAP_lens),
             D_Wf0=data$d_vectorspp,
             D_LensA=data$d_AAP[1,], # units: hours
             D_LensL=data$d_LAP[1,], # units: hours
             D_LensI=data$d_IAP[1,], # units: hours
             D_RepsA=data$d_AAP[2,],
             D_RepsL=data$d_LAP[2,],
             D_RepsI=data$d_IAP[2,],
             D_InfsA=data$d_AAP[3,],
             D_InfsL=data$d_LAP[3,],
             D_InfsI=data$d_IAP[3,],
             D_bgLens=data$d_durations, # units: hours
             D_lsPars=c(D_LSin,D_numPtsPdin))

  # Note: options() both sets and returns options...
  # so that when one calls eg options(mc.cores = 4), it both sets mc.cores = 4 AND returns list of previous values.
  old_opts <- options(mc.cores = if (mc.parallel == 1) parallel::detectCores() else 1)   # Set mc.cores based on mc.parallel argument
  on.exit(options(old_opts), add = TRUE)

  #### run stan on the data as a model OF AP assays for PT data ################
  #model=readRDS(system.file("stan", "APmodel_PT_virus.rds", package = "EpiPv"))

  fitB1 = rstan::sampling(stanmodels$APmodel_PT_virus, data = c(dat3,D_ppc=0),
                          iter = numIter,
                          warmup = numWarm,
                          chains=numChains,
                          control = list(max_treedepth = treeDepth,stepsize=0.01,adapt_delta=adaptVal))

  rstan::check_energy(fitB1)

  draws_array <- posterior::as_draws_array(fitB1, inc_warmup = FALSE)
  params<- posterior::subset_draws(draws_array, variable = c("al[1]", "be[1]", "mu_lat[1]", "lat[1]", "bD[1]","lp__"))
  parDim=dim(params)

  params1<- posterior::subset_draws(draws_array, variable = c("al[1]", "be[1]", "mu[1]", "lat[1]"))

  #### Compute the diagnostics from the draws object
  summary_stats <- posterior::summarise_draws(draws_array)

  #### Retrieve the convergence results for extra safety
  sampler_params <- rstan::get_sampler_params(fitB1, inc_warmup = FALSE)

  dts=numeric(0)
  trds=numeric(0)
  for (kkk in 1:numChains) {
    dts=c(dts,sum(sampler_params[[kkk]][, "divergent__"]))
    trds=c(trds,max(sampler_params[[kkk]][, "treedepth__"]))
  }

  diagntcs <- list(
    divergent_transitions = sum(dts),
    max_treedepth_exceeded = (max(trds)>= treeDepth)
  )


  message(
    "MESSAGE: Check core convergence diagnostics and model fit!\n",
    "- First: Review RStan sampling messages for warnings (do NOT suppress them).\n",
    "- R-hat < 1.05 for all parameters - see array 2.\n",
    "- Effective Sample Size (ESS) ~>= 100 per chain - see array 2.\n",
    "- Posterior predictive fit: compare simulated vs. observed (arrays 3, 4, 5).\n",
    "- No divergent transitions or treedepth exceedances - see converge_results.\n",
    "More details:\n",
    "   Function help & vignette (includes advanced diagnostics)\n",
    "   Stan diagnostics guide: https://mc-stan.org/learn-stan/diagnostics-warnings.html\n"
  )


  # for validation plots
  # gather credible intervals from stan generated quantities (corresponds to reproduction of exp data)
  simA=array(rep(0,2*dim(data$d_AAP)[2]),dim=c(dim(data$d_AAP)[2],2))
  simL=array(rep(0,2*dim(data$d_LAP)[2]),dim=c(dim(data$d_LAP)[2],2))
  simI=array(rep(0,2*dim(data$d_IAP)[2]),dim=c(dim(data$d_IAP)[2],2))

  y_sim_A=array(rep(0,parDim[1]*parDim[2]*dim(data$d_AAP)[2]),dim=c(parDim[c(1,2)],dim(data$d_AAP)[2]))
  y_sim_L=array(rep(0,parDim[1]*parDim[2]*dim(data$d_LAP)[2]),dim=c(parDim[c(1,2)],dim(data$d_LAP)[2]))
  y_sim_I=array(rep(0,parDim[1]*parDim[2]*dim(data$d_IAP)[2]),dim=c(parDim[c(1,2)],dim(data$d_IAP)[2]))

  for (sss in 1:(dim(data$d_AAP)[2])) {
    simA[sss,]=rstan::summary(fitB1,paste0('y_simul_A[',sss,']'))$summary[c(4,8)]
    y_sim_A[,,sss]<- posterior::subset_draws(draws_array, variable = paste0('y_simul_A[',sss,']'))
  }
  for (ggg in 1:(dim(data$d_LAP)[2])) {
    simL[ggg,]=rstan::summary(fitB1,paste0('y_simul_L[',ggg,']'))$summary[c(4,8)]
    y_sim_L[,,ggg]<- posterior::subset_draws(draws_array, variable = paste0('y_simul_L[',ggg,']'))
  }
  for (mmm in 1:(dim(data$d_IAP)[2])) {
    simI[mmm,]=rstan::summary(fitB1,paste0('y_simul_I[',mmm,']'))$summary[c(4,8)]
    y_sim_I[,,mmm]<- posterior::subset_draws(draws_array, variable = paste0('y_simul_I[',mmm,']'))
  }

  #### for Gelman Bayesian R2
  var_simA <- apply(y_sim_A, c(1,2), var)  # Variance across iterations and chains for each observation (assay_length)
  var_simL <- apply(y_sim_L, c(1,2), var)  # Variance across observations for each MCMC sample
  var_simI <- apply(y_sim_I, c(1,2), var)  #

  # Compute residuals (y_rep - y), where simA and simI are 3D
  # Since dat4b$D_InfsA and dat4b$D_InfsI are vectors (1 x assay_length), repeat them over chains and iterations
  residualsA <- y_sim_A - array(rep(dat3$D_InfsA, each = (numChains * (numIter-numWarm))), dim = dim(y_sim_A))
  residualsL <- y_sim_L - array(rep(dat3$D_InfsL, each = (numChains * (numIter-numWarm))), dim = dim(y_sim_L))
  residualsI <- y_sim_I - array(rep(dat3$D_InfsI, each = (numChains * (numIter-numWarm))), dim = dim(y_sim_I))

  # Compute variance of residuals across iterations and chains
  var_residualA <- apply(residualsA, c(1,2), var)
  var_residualL <- apply(residualsL, c(1,2), var)
  var_residualI <- apply(residualsI, c(1,2), var)

  # Compute Bayesian R^2
  r2_bayes_A <- var_simA / (var_simA + var_residualA)
  r2_bayes_L <- var_simL / (var_simL + var_residualL)
  r2_bayes_I <- var_simI / (var_simI + var_residualI)

  df1=list(lenA=dat3$D_LensA,
          propA=dat3$D_InfsA/dat3$D_RepsA,
          simulLA=simA[,1]/dat3$D_RepsA,
          simulUA=simA[,2]/dat3$D_RepsA)

  df2=list(lenL=dat3$D_LensL,
          propL=dat3$D_InfsL/dat3$D_RepsL,
          simulLL=simL[,1]/dat3$D_RepsL,
          simulUL=simL[,2]/dat3$D_RepsL)

  df3=list(lenI=dat3$D_LensI,
          propI=dat3$D_InfsI/dat3$D_RepsI,
          simulLI=simI[,1]/dat3$D_RepsI,
          simulUI=simI[,2]/dat3$D_RepsI)

  df4=list(bayesR2_mn=c(mean(r2_bayes_A),mean(r2_bayes_L),mean(r2_bayes_I)),
           bayesR2_sd=c(sd(r2_bayes_A),sd(r2_bayes_L),sd(r2_bayes_I)))


  return(list(array0 = params, array1 = params1, array2 = summary_stats[,variable = c("rhat", "ess_bulk, ","ess_tail")], array3=df1, array4=df2, array5=df3, array6=df4,converge_results=diagntcs));

}
