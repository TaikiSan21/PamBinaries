#' @title Read Whistle & Moan Detector Header
#' 
#' @description Reads file header information specific to the 
#'   Whistle & Moan Detector module
#'   
#' @param file binary file to be read
#' 
#' @return header information for the Whistle & Moan Detector module
#' 
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#' 
readWMDHeader <- function(file) {
    header <- readStdModuleHeader(file)
    if((header$binaryLength != 0) & (header$version >= 1)) {
        header$delayScale <- pamBinRead(file, 'int32', n=1)
    }
    return(header)
}