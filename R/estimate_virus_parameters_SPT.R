#' Estimate SPT virus parameters using Bayesian inference with rstan
#'
#' Runs MCMC sampling using a precompiled Stan model for SPT plant virus. Analyses SPT access period data directly estimating 4 parameters:
#' (\code{c1[1]}, \code{c3[1]}, \code{c2[1]},\code{bD[1]}, the first 3 are composite and the 4th corresponds to insect lab mortality rate).
#' The composite parameters can be unpacked into \code{al[1]}, \code{be[1]}, \code{mu[1]} - acquisition, inoculation, and insect recovery rates.
#'
#' @param data A list containing the AP experiment data for Stan (required).
#' @param D_LSin Upper bound on insect vector lifespan in days, sets the vector survival discretisation (optional, default = 50).
#' @param D_numPtsPdin Number of points per day of insect vector lifespan, sets the vector survival discretisation (optional, default = 1).
#' @param mcmcOptions A numeric vector of length 2:
#'   The first element specifies the number of warm-up iterations (optional, default = 500),
#'   and the second element specifies the number of post-warm-up iterations (optional, default = 1000).
#' @param numChainsIn Numeric: number of Markov chains (optional, default = 4).
#' @param mc.parallel Binary: whether to use parallelisation for sampling (optional, default = 0, i.e. 1 core only).
#' @return A list with seven elements:
#'   \describe{
#'     \item{array0}{First set of MCMC chains for estimated parameters (\code{c1[1]}, \code{c3[1]}, \code{c2[1]},\code{bD[1]} - and additionally \code{lp__} for reference only).}
#'     \item{array1}{Second set of MCMC chains for virus parameters (\code{al[1]}, \code{be[1]}, \code{mu[1]}).}
#'     \item{array2}{summary_stats (rhat, ess_bulk, ess_tail). \code{[1]}: R-hat statistic for convergence (should be close to 1); \code{[2-3]}: statistics for n-eff.}
#'     \item{array3}{Validation list, AAP input test plant infection data, forward simulated 2.5th percentile, forward simulated 97.5th percentile.}
#'     \item{array4}{Validation list, IAP input test plant infection data, forward simulated 2.5th percentile, forward simulated 97.5th percentile.}
#'     \item{array5}{Mean and sd Bayesian R-squared values for model fit assessment, for AAP and IAP assays.}
#'     \item{converge_results}{A list containing key sampler diagnostics: divergent transitions, maximum tree depth exceedances,
#'                      See also screen print for acceptability of energy Bayesian fraction of missing information (E-BFMI).).}
#'   }
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
#' 4. **Posterior predictive fit (simulated_data)** - Compare model-simulated data with the observed data:
#'    - \code{array3}: AAP forward simulation 95% credible intervals and original data
#'    - \code{array4}: IAP forward simulation 95% credible intervals and original data
#'    (\href{https://mc-stan.org/docs/stan-users-guide/posterior-predictive-checks.html}{Stan posterior predictive checks}).
#'
#' **Additional diagnostics for troubleshooting**
#'
#' These are useful when core checks indicate problems, or for deeper exploration of model behaviour:
#'
#' - **Bayesian R^2** (\code{r2_bayesA}, \code{r2_bayesI}) - Measures explanatory power;
#'   low values suggest poor fit to the data (i.e., 2 values for AAP assay, IAP assay).
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
#' data("ap_data_sim_SPT", package = "EpiPvr")
#' # run low warm-up and iterations (mcmcOptions) for quick example only
#' EVSPT_pub=estimate_virus_parameters_SPT(
#'   data=ap_data_sim_SPT,D_LSin=5,D_numPtsPdin=1,
#'   mcmcOptions=c(25,50),numChainsIn=1,mc.parallel=0
#' )
#' @export
estimate_virus_parameters_SPT <- function(data,D_LSin=50,D_numPtsPdin=1,mcmcOptions=c(500,1000),numChainsIn=4,mc.parallel=0){


  numWarm=mcmcOptions[1]
  numIter=mcmcOptions[2]
  adaptVal=0.95 #1-(10^-1)
  treeDepth=10 #15
  numChains=numChainsIn

  message(paste('numChains =',numChainsIn))

  ##### diags should be negative in assay durations matrix (corresponds to separate duration vector that is varied) #####
  if (sum(diag(data$d_durations)>0)!=0){
    stop('input error! see data d_durations');
  }

  ##### data must have acquisitiona access and inoculationa access to be valid AP data #####
  if (any((data$d_durations[1,2]==0),(data$d_durations[2,1]==0))){
    stop('input error! see data d_durations');
  }

  ##### data should only involve SPT virus #####
  if (!identical(data$d_virusType,"SPT")){
    stop('input error! see data d_virusType');
  }

  # assembling all the data inputs for the stan estimation
  dat4b= list(D_NumGrps=c(dim(data$d_AAP)[2],dim(data$d_IAP)[2]),#length(AAP_lens),length(LAP_Reps_mg)
              D_Wf0=data$d_vectorspp,
              D_LensA=data$d_AAP[1,], # units: hours
              D_LensI=data$d_IAP[1,], # units: hours
              D_RepsA=data$d_AAP[2,],
              D_RepsI=data$d_IAP[2,],
              D_InfsA=data$d_AAP[3,],
              D_InfsI=data$d_IAP[3,],
              D_bgLens=data$d_durations, # units: hours
              D_lsPars=c(D_LSin,D_numPtsPdin))

  # Note: options() both sets and returns options...
  # so that when one calls eg options(mc.cores = 4), it both sets mc.cores = 4 AND returns list of previous values.
  old_opts <- options(mc.cores = if (mc.parallel == 1) parallel::detectCores() else 1)   # Set mc.cores based on mc.parallel argument
  on.exit(options(old_opts), add = TRUE)

  #### run stan on the data as a model OF AP assays for SPT data ################
  #model=readRDS(system.file("stan", "APmodel_SPT_virus.rds", package = "EpiPv"))


  fitB2 = rstan::sampling(stanmodels$APmodel_SPT_virus, data = dat4b,
                          iter = numIter,
                          warmup = numWarm,
                          chains=numChains,
                          control = list(max_treedepth = treeDepth,adapt_delta=adaptVal,stepsize=0.01))

  rstan::check_energy(fitB2)

  #### return chains for main parameters ####
  draws_array <- posterior::as_draws_array(fitB2, inc_warmup = FALSE)

  params<- posterior::subset_draws(draws_array, variable = c("c1[1]","c3[1]","c2[1]","bD[1]","lp__"))

  params1<- posterior::subset_draws(draws_array, variable = c("al[1]", "be[1]", "mu[1]"))

  parDim=dim(params)

  #### Compute the diagnostics from the draws object
  summary_stats <- posterior::summarise_draws(draws_array)

  #### Retrieve the convergence results for extra safety
  sampler_params <- rstan::get_sampler_params(fitB2, inc_warmup = FALSE)

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
    "- Posterior predictive fit: compare simulated vs. observed (arrays 3 & 4).\n",
    "- No divergent transitions or treedepth exceedances - see converge_results.\n",
    "More details:\n",
    "   Function help & vignette (includes advanced diagnostics)\n",
    "   Stan diagnostics guide: https://mc-stan.org/learn-stan/diagnostics-warnings.html\n"
  )

  #### for validation plots
  # gather credible intervals from stan generated quantities (corresponds to reproduction of exp data)
  simA=array(rep(0,2*dim(data$d_AAP)[2]),dim=c(dim(data$d_AAP)[2],2))
  simI=array(rep(0,2*dim(data$d_IAP)[2]),dim=c(dim(data$d_IAP)[2],2))

  y_sim_A=array(rep(0,parDim[1]*parDim[2]*dim(data$d_AAP)[2]),dim=c(parDim[c(1,2)],dim(data$d_AAP)[2]))
  y_sim_I=array(rep(0,parDim[1]*parDim[2]*dim(data$d_IAP)[2]),dim=c(parDim[c(1,2)],dim(data$d_IAP)[2]))

  for (sss in 1:(dim(data$d_AAP)[2])) {
    simA[sss,]=rstan::summary(fitB2,paste0('y_simul_A[',sss,']'))$summary[c(4,8)]
    y_sim_A[,,sss]=posterior::subset_draws(draws_array, variable = paste0('y_simul_A[',sss,']'))
  }
  for (mmm in 1:(dim(data$d_IAP)[2])) {
    simI[mmm,]=rstan::summary(fitB2,paste0('y_simul_I[',mmm,']'))$summary[c(4,8)]
    y_sim_I[,,mmm]<- posterior::subset_draws(draws_array, variable = paste0('y_simul_I[',mmm,']'))
  }

  # Compute variance across iterations and chains for each observation (across the assay length)
  var_simA <- apply(y_sim_A, c(1,2), var)  # Variance across iterations and chains for each observation (assay_length)
  var_simI <- apply(y_sim_I, c(1,2), var)  # Same for portion B

  # Compute residuals (y_rep - y), where simA and simI are 3D
  # Since dat4b$D_InfsA and dat4b$D_InfsI are vectors (1 x assay_length), repeat them over chains and iterations
  residualsA <- y_sim_A - array(rep(dat4b$D_InfsA, each = (numChains * (numIter-numWarm))), dim = dim(y_sim_A))
  residualsI <- y_sim_I - array(rep(dat4b$D_InfsI, each = (numChains * (numIter-numWarm))), dim = dim(y_sim_I))

  # Compute variance of residuals across iterations and chains
  var_residualA <- apply(residualsA, c(1,2), var)
  var_residualI <- apply(residualsI, c(1,2), var)

  # Compute Bayesian R^2
  r2_bayes_A <- var_simA / (var_simA + var_residualA)
  r2_bayes_I <- var_simI / (var_simI + var_residualI)

  ### collating data for output objects in list
  df1=list(lenA=dat4b$D_LensA,
          propA=dat4b$D_InfsA/dat4b$D_RepsA,
          simulLA=simA[,1]/dat4b$D_RepsA,
          simulUA=simA[,2]/dat4b$D_RepsA)

  df2=list(lenI=dat4b$D_LensI,
           propI=dat4b$D_InfsI/dat4b$D_RepsI,
           simulLI=simI[,1]/dat4b$D_RepsI,
           simulUI=simI[,2]/dat4b$D_RepsI)

  df3=list(bayesR2_mn=c(mean(r2_bayes_A),mean(r2_bayes_I)),
           bayesR2_sd=c(sd(r2_bayes_A),sd(r2_bayes_I)))

  return(list(array0 = params, array1 = params1, array2 = summary_stats[,variable = c("rhat", "ess_bulk, ","ess_tail")], array3=df1, array4=df2, array5=df3, converge_results=diagntcs));
}
