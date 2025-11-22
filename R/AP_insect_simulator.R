#' Simulate AP assay data for individual insects
#'
#' Simulates AP experiment data by modeling each insect's virus acquisition, latency, recovery, and inoculation independently.
#'
#' @param lmark_in Numeric vector of length 3 specifying assay durations.
#' @param smarkpams_in Numeric vector of length 4 specifying virus parameters (rates per hour):
#'   The first element specifies the acquisition rate,
#'   the second the inoculation rate,
#'   the third the latent progression rate; code -1 for SPT
#'   the forth the insect recovery rate
#' @param isVerbose Logical (default: FALSE). If TRUE, prints information from nested functions.
#' @param virusType Character. Plant virus type, either SPT or PT.
#' @return A numeric vector of length 2:
#'   \describe{
#'     \item{\code{[1]}}{Binary infection status of the test plant (1 = infected, 0 = not infected).}
#'     \item{\code{[2]}}{Number of virus re-acquisition events (relevant for SPT viruses).}
#'   }
#' @keywords internal
#' @details
#' AP insect simulator models insects in AP assays as acquiring virus, progressing through latency, recovering, and inoculating independently.
#' Since acquisition can happen in the initial period, latency progression may occur either within the acquisition period, during the latency-specific period,
#' or in the inoculation period, affecting the final inoculation time. The function outputs whether an insect successfully inoculates a test plant
#' and the number of virus re-acquisition events (for SPT viruses).
AP_insect_simulator <- function(lmark_in,smarkpams_in,isVerbose,virusType) {

  alrate=smarkpams_in[1]
  berate=smarkpams_in[2]
  gamrate=smarkpams_in[3]
  murate=smarkpams_in[4]

  lmark=lmark_in

  if (isVerbose){
    message(paste('read-in assay lengths: ',lmark))
    message(paste('read-in rATES: ',smarkpams_in))
  }

  ##### data for SPT virus not expected to have latent progression rate #####
  if (identical(virusType,"SPT")&&(!is.na(gamrate))){
    stop('input error! SPT virus should have NA gamma');
  }

  # SIMULATE TIMES UNTIL ACQUISITION, LATENCY PROGRESSION, RECOVERY
  tA_sim=gen_exp(1,alrate,virusType)
  tR_sim=gen_exp(1,murate,virusType)
  if (is.na(gamrate)) {
    tL_sim=gen_exp(1,Inf,virusType)
  }else{
    tL_sim=gen_exp(1,gamrate,virusType)
  }

  # provide new acquisition, latency progression, and recovery times if there is recovery before the end of the AAP
  numReacq=0
  while((tA_sim+tL_sim+tR_sim)<lmark[1]){      #resample
    tA_sim=gen_exp(1,alrate,virusType)+(tA_sim+tL_sim+tR_sim)
    tR_sim=gen_exp(1,murate,virusType)
    if (is.na(gamrate)) {
      tL_sim=gen_exp(1,Inf,virusType)
    }else{
      tL_sim=gen_exp(1,gamrate,virusType)
    }
    numReacq=numReacq+1
    if (isVerbose){
      message(paste('...Re-acquisition...'))
    }
  }
  tmp_Durations=inoc_durtn_calculator(lmark,c(tA_sim,tL_sim,tR_sim),isVerbose)

  # ENACT inoculation
  rnum=gen_exp(1,berate,virusType)
  TPI=ifelse((rnum<tmp_Durations),1,0)
  return(c(TPI,numReacq))
}


#' Safely generate random exponential data
#'
#' A helper function only needed within AP_insect_simulator that generates random numbers from an exponential distribution
#' using the base/stats R function `stats::runif`.
#' This function simulates a set of random numbers with a fixed rate safely (ie for different input cases stemming from virus types).
#'
#' @param n Numeric: number of points to sample.
#' @param rin Numeric: event rate.
#' @param virusType Character: SPT or PT plant virus type.
#' @return Numeric vector: n random exponential number.
#' @keywords internal
gen_exp <- function(n, rin, virusType) {
  if (is.nan(rin)) {
    stop('undefined exponential rate')
  } else if(is.na(rin)){
    stop('missing exponential rate')
  }else if(rin==Inf){
    out=rep(0,n)
    # GENERATE a warning and test this warning! but only for PT - inf rin is expected for SPT
    if (virusType=="SPT") {
      # n.b 'rin==Inf' is an acceptable case for SPT (to cover progression through latency which is assumed instant for SPT virus)
    }else{
      warning('infinite exponential was specified - proceeding with zero time sample', call. = TRUE, immediate. = TRUE, noBreaks. = FALSE,
              domain = NULL)
    }
  }
  else if(rin==0){
    out=rep(Inf,n)
    warning('exponential rate of zero was specified - proceeding with infinite time sample', call. = TRUE, immediate. = TRUE, noBreaks. = FALSE,
            domain = NULL)
  }else{
    if (rin <= 0){
      stop("Rate parameter must be positive")
    }else{
      out=-log(stats::runif(n))/rin  # Inverse sample from cdf
    }
  }
  return(out)
}


