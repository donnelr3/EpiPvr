#'
#'  Helper function for effective inoculation period calculator
#'
#' @param lmark_in Numeric vector of length 3 specifying assay durations.
#' @param sim_Ts_in Numeric vector of length 4 specifying sampled event times (hours) related to virus parameters (rates hr^-1) (rates per hour)
#' @param isitVerbose Logical (default: `FALSE`). If `TRUE`, prints information from nested functions.
#' @return Numeric: effective inoculation duration.
#' @keywords internal
inoc_durtn_calculator <- function(lmark_in,sim_Ts_in,isitVerbose) {
  
  tA_sim=sim_Ts_in[1] # little t_'s are event times
  tL_sim=sim_Ts_in[2]+tA_sim
  tR_sim=sim_Ts_in[3]+tL_sim
  T_Ain=lmark_in[1] # big T_'s are assay timepoints
  T_Lin=lmark_in[2]+T_Ain
  T_Iin=lmark_in[3]+T_Lin
  if (isitVerbose){
    if (T_Ain<tA_sim){message('warning no acquisition in acquisition period')}
    if (T_Lin<tL_sim){message('warning no latency pass in latent period but can still pass in inoc period')}
    if (T_Lin>tR_sim){message('warning loss of pathogen prior to inoculation period')}
  }
  if ((tA_sim>T_Ain)||((tL_sim)>T_Iin)||((tR_sim)<T_Lin)) {
    output=0
  }else{
    output=max((T_Iin-T_Lin)-max(0,(T_Iin-(tR_sim)))-max(0,(tL_sim-T_Lin)),0)
  }
  return(output)
}
