#' @title Read Click Trigger Header
#' 
#' @description Reads file header information specific to the click trigger module
#'   
#' @param file binary file to be read
#' 
#' @return header information for the click trigger
#' 
#' @author Taiki Sakai \email{taiki.sakai@noaa.gov}
#' 
readClickTriggerHeader <- function(file) {
    header <- readStdModuleHeader(file)
    if(header$binaryLength != 0) {
        header$channelMap <- pamBinRead(file, 'int32', n=1)
        nChan <- countChannels(header$channelMap)
        header$calibration <- pamBinRead(file, 'float', n=nChan)
    } else {
        header$calibration <- NULL
    }
    return(header)
}