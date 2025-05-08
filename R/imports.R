#' @importFrom sp SpatialPoints spDists CRS
#' @importFrom spacetime STIDF
#' @importFrom trajectories Track frechetDist
#' @importFrom purrr map_df map_dfr map
#' @importFrom magic adiag
#' @importFrom utils combn globalVariables
#' @importFrom stats quantile sd qnorm
#' @importFrom dplyr %>%
#' @importFrom future plan
#' @importFrom parallelly availableCores
#' @importFrom progressr progressor with_progress handlers
#' @importFrom stats pnorm
#' @import parallel dplyr
utils::globalVariables(c("ID","trip","ID1b", "trip1", "ID2b","trip2","ID1","ID2","d","ID2b","trip2b","i","j","idtrip"))
NULL
