#' Main function for AP assay simulation
#'
#' This function simulates the assay data of an AP experiment (calling the helper function AP_insect_simulator.r to simulate on a per insect basis).
#'
#' @param assayinput_in Numeric array with 3 rows (required): \code{[1,]} assay variable durations, \code{[2,]} reps per duration, \code{[3,]} zero entries for infected test plants (to be populated with simulated values).
#' @param fixedComponent_in Numeric vector of 3 fixed durations in the order AAP, LAP, IAP; nb. one of these 3 will not be relevant for the sub-assay as it is varied in the assay (code -1) e.g. if AAP sub-assay \code{[1]} is -1, if LAP sub-assay \code{[2]} is -1, if IAP sub-assay \code{[3]} is -1 (required).
#' @param WF_in Number of insect vectors in cohort (required).
#' @param smarkpams_in Numeric vector: 4 elements representing virus rate parameters (hr^-1) (required):
#'   \code{[1]}, acquisition rate;
#'   \code{[2]}, inoculation rate;
#'   \code{[3]}, latent progression rate WITH code -1 for SPT (no significant latency for SPT virus);
#'   \code{[4]}, insect recovery rate.
#' @param isVerbose Binary (optional, default = 0). If 1, prints information from nested functions.
#' @param virusTypeIn Character. Plant virus type, either `"SPT"` or `"PT"`  (required).
#' @return A numeric array with 3 rows corresponding to assay structure and outcome:
#'   The first row specifies the assay variable durations,
#'   the second the reps per duration,
#'   the third the simulated number of infected test plants.
#'
#' @examples
#' assay1=AP_assay_simulator(
#'   assayinput_in=rbind(c(6,7,8),rep(10,3),rep(0,3)),
#'   fixedComponent_in=c(10,-1,10),WF_in=5,
#'   smarkpams_in=rep(0.1,4),isVerbose=0,virusTypeIn='PT'
#' )
#' @export
AP_assay_simulator <- function(assayinput_in, fixedComponent_in, WF_in, smarkpams_in,isVerbose=0,virusTypeIn) {

  if (any(is.na(assayinput_in))||any(is.nan(assayinput_in)))
    stop('assay input contains missing or NaN values')
  else if (any(is.na(fixedComponent_in))||any(is.nan(fixedComponent_in)))
    stop('assay durations contains missing or NaN values')
  else if (any(is.na(WF_in))||any(is.nan(WF_in)))
    stop('num insects contains missing or NaN values')

  fixedComponent=fixedComponent_in
  assay_dim=dim(assayinput_in)
  T_vec=assayinput_in[1,]
  R_vec=assayinput_in[2,]
  I_vec=assayinput_in[3,]

  ### ASSAY
  for (nv in 1:length(T_vec)) {
    inoc_byRepsbyWf=matrix(0,WF_in,R_vec[nv])
    for (nwr in 1:R_vec[nv]) {
      cntrBig=0
      fixedComponent[fixedComponent_in==-1]=T_vec[nv] # '-1' is denoting the varying part of the assay. Here replacing with each instance
      for (nwf in 1:WF_in) {
        tmp=AP_insect_simulator(fixedComponent,smarkpams_in,isVerbose,virusTypeIn)
        inoc_byRepsbyWf[nwf,nwr]=tmp[1]               # filling up a binary test plant inoculation matrix for each rep and each insect
        cntrBig=cntrBig+tmp[2]                        # keeping track of number of re-acquisitions (requires recovery to occur within AAP) mainly for testing
      }
    }
    I_vec[nv]=sum((colSums(inoc_byRepsbyWf)>0)*1)        #  checking that at least one of the insects inoculated (i.e., test plant infection)
  }
  return(rbind(T_vec,R_vec,I_vec))
}
