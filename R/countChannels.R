#' @title Count Number of Active Channels
#' 
#' @description Counts the number of active channels
#'   given a channel mapping
#'   
#' @param channelMap Mapping of channels as a binary number
#' 
#' @return The number of active channels (number of ones)
#' 
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#' 
#' @note Altered from original script to loop through 30 instead
#'   32 because R stores only 32 bit integers. Should not ever have
#'   enough channels for this to matter.
#'   
countChannels <- function(channelMap) {
    nC <- 0
    j <- 1
    for(i in 1:30) {
        if(bitwAnd(channelMap, j) != 0) {
            nC <- nC + 1
        }
        j <- j * 2
    }
    return(nC)
}
