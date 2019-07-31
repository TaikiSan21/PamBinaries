#' @title Read LTSA Header
#' 
#' @description Reads file header information specific to the LTSA module
#'   
#' @param file binary file to be read
#' 
#' @return header information for the LTSA module
#' 
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#' 
readLTSAHeader <- function(file) {
    header <- readStdModuleHeader(file)
    if(header$binaryLength != 0) {
        header$fftLength <- pamBinRead(file, 'int32', n=1)
        header$fftHop <- pamBinRead(file, 'int32', n=1)
        header$intervalSeconds <- pamBinRead(file, 'int32', n=1)
    }
    return(header)
}