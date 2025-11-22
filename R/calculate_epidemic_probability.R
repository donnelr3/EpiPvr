#' Main function for calculating epidemic probability
#'
#' This function calculates epidemic probability from viral and local parameters according to inoculum state (calling the helper function solveInoculumStatesBP.r to calculate transition and fertility matrices)
#'
#' @param numberInsects Numeric: number of insect vectors per plant (required).
#' @param start_intervl Numeric: initial guess for expected time in hours until next event (optional, default = 0.1). This is fine-tuned from within the function. Needed to calculate efficient time step i.e. largest value such that sum of probabilities of individual events < 1.
#' @param localParameters A numeric vector of length 5 corresponding to event rates (day^-1) (required):
#'   The first element specifies the Vector dispersal rate,
#'   the second the Roguing rate,
#'   the third the Vector mortality rate,
#'   the forth the Harvesting rate,
#'   the fifth the Plant latent progression rate
#' @param virusParameters Numeric vector of length 3 specifying virus parameters (day^-1) (required):
#'   \describe{
#'     \item{\code{[1]}}{Acquisition rate.}
#'     \item{\code{[2]}}{Inoculation rate.}
#'     \item{\code{[3]}}{Insect recovery rate.}
#'   }
#' @param thresh Numeric: equilibrium condition for extinction probability (optional, default = 10^(-7)). If no inoculum state changes by more than thresh in a timestep then function assumes extinction probability has stabilised.
#' @return Numeric vector of length n: epidemic probability for each of n inoculum states (where n=(numberInsects+1)*3-1).
#'
#' **Preprint reference**
#'
#' The models implemented in this function follow Donnelly et al. (2025, preprint), originally implemented in the EpiPv GitHub package.
#'
#' @references
#' Donnelly, R., Tankam, I. & Gilligan, C. (2025).
#' "Plant pathogen profiling with the EpiPv package."
#' EcoEvoRxiv, \doi{10.32942/X29K9P}.
#'
#' When available, please cite the published version.
#'
#' @examples
#' qm_out=calculate_epidemic_probability(
#'     numberInsects=3,start_intervl=(1/10),
#'     localParameters=rep(0.1,5),
#'     virusParameters=rep(0.1,3),thresh=10^(-7)
#' )
#' @export
calculate_epidemic_probability <- function(numberInsects,start_intervl=(1/10),localParameters,virusParameters,thresh=10^(-7)){

  numVariables <- ((numberInsects+1)*3)-1; # inoculum states determined by # insects per plant
  #### INITIAL section selects an efficient timestep (because the prob of ANY event in a time step must be <=1... hence must calculate optimal time step)
  interval_indN=1/start_intervl # total rate, increases by 1 in search loop; as we search for the smallest time step such that prob no event occurs is non-zero
  interval_rec=24*interval_indN
  # solveInoculumStatesBP computes the transition matrix between inoculum states -  initially call this repeatedly to figure out optimal time step
  result <- solveInoculumStatesBP(numVariables, localParameters, virusParameters, 1/interval_rec)
  dr1=dim(result)[1]
  ttilde=result[1:(numVariables+1),]
  f=result[(numVariables+1+1):dr1,]
  # the best choice of interval length is the biggest which produces only positive entries in ttilde - see MS appendix S3
  while (sum(ttilde<0)) {
    interval_indN=interval_indN+1
    interval_rec=24*interval_indN
    result=solveInoculumStatesBP(numVariables, localParameters, virusParameters, 1/interval_rec)
    dr1=dim(result)[1]
    ttilde=result[1:(numVariables+1),] # transition matrix
    f=result[(numVariables+1+1):dr1,] # fertility matrix
  }
  #message(paste0('efficient interval length in mins is ',24*60*(1/interval_rec)))

  #### SUBSEQUENT section iterates branching process method for finding extinction probability, see ms appendix S3
  itimeCtr <- 0
  k <- dim(ttilde)[2]
  qm <- rep(0, k)           # extinction probability vector
  qprime <- rep(0, k)      # updated extinction probability
  qm_old=rep(-1, k)

  # Iterating time loop
  #message('start extinction solve')
  while(sum(abs(qm-qm_old)>thresh)>0){#(itimeCtr<100000){#
    qm_old=qm
    itimeCtr=itimeCtr+1
    gb <- 1 - f[numVariables-numberInsects+1, ] + f[numVariables-numberInsects+1, ] * qm[numVariables-numberInsects+1]  # Bernoulli term
    gt <- t(ttilde) %*% c(qm, 1)     # Matrix multiplication (nx(n+1)) %*% ((n+1)x1)
    qprime <- gt * gb  # update extinction probability
    qm <- qprime  # Update q with qprime for next iteration
  }


  #message('end extinction solve')
  return(1-qm); # since branching procseses calculation is based upon extinction probability, convert to epidemic probability epi prob=1-ext prob
}
