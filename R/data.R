#' Simulated access period experiment Data for a hypothetical SPT virus
#'
#' This dataset contains simulated AP assay results for an SPT virus (AAP and IAP sub-assays)
#'
#' @format A List that is a collection of 5 data objects and 3 attributes
#' \describe{
#'   \item{d_AAP}{Numeric 3 x n array for acquistion access sub-assay: variable duration; reps per duration, no. infected test plants per duration}
#'   \item{d_IAP}{Numeric 3 x n array for inoculation access sub-assay: variable duration; reps per duration, no. infected test plants per duration}
#'   \item{d_durations}{Numeric 2 x 2 assay of fixed durations (FD): row 1, FDs for acquistion access sub-assay; row 2, FDs acquistion access sub-assay. nb  diagonal entries must be '-1' indicating redundancy since these elements correspond to the variable durations}
#'   \item{d_vectorspp}{Numeric number of insects per plant}
#'   \item{d_pathogen}{Character representing underlying virus}
#' }
#' #' \describe{
#'   \item{"alpha"}{Numeric: acquisition rate (h^-1); used to generate the data or '-1' if unknown}
#'   \item{"beta"}{Numeric: inoculation rate (h^-1); used to generate the data or '-1' if unknown}
#'   \item{"mu"}{Numeric: insect recovery rate (h^-1); used to generate the data or '-1' if unknown}
#' }
#' @source Simulated using `AP_assay_simulator()`
"ap_data_sim_SPT"

#' Simulated access period experiment Data for a hypothetical PT virus
#'
#' This dataset contains simulated AP assay results for a PT virus (AAP LAP and IAP sub-assays)
#'
#' @format A List that is a collection of 6 data objects and 4 attributes
#' \describe{
#'   \item{d_AAP}{Numeric 3 x n array for acquistion access sub-assay: variable duration; reps per duration, no. infected test plants per duration}
#'   \item{d_LAP}{Numeric 3 x n array for latent access sub-assay: variable duration; reps per duration, no. infected test plants per duration}
#'   \item{d_IAP}{Numeric 3 x n array for inoculation access sub-assay: variable duration; reps per duration, no. infected test plants per duration}
#'   \item{d_durations}{Numeric 3 x 3 assay of fixed durations (FD): row 1, FDs for acquistion access sub-assay; row 2, FDs latent access sub-assay; row 3, FDs inoculation access sub-assay. nb  diagonal entries must be '-1' indicating redundancy since these elements correspond to the variable durations}
#'   \item{d_vectorspp}{Numeric number of insects per plant}
#'   \item{d_pathogen}{Character representing underlying virus}
#' }
#' #' \describe{
#'   \item{"alpha"}{Numeric: acquisition rate (h^-1); used to generate the data or '-1' if unknown}
#'   \item{"beta"}{Numeric: inoculation rate (h^-1); used to generate the data or '-1' if unknown}
#'   \item{"gamma"}{Numeric: latent progression rate (h^-1); used to generate the data or '-1' if unknown}
#'   \item{"mu"}{Numeric: insect recovery rate (h^-1); used to generate the data or '-1' if unknown}
#' }
#' @source Simulated using `AP_assay_simulator()`
"ap_data_sim_PT"