#' @title Read Noise Monitor Header
#' 
#' @description Reads file header information specific to the Noise Monitor module
#'   
#' @param file binary file to be read
#' 
#' @return header information for the Noise Monitor module
#' 
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#' 
readNoiseMonHeader <- function(file) {
    header <- readStdModuleHeader(file)
    if(header$binaryLength != 0) {
        header$nBands <- pamBinRead(file, 'int16', n=1)
        header$statsTypes <- pamBinRead(file, 'int16', n=1)
        header$loEdges <- pamBinRead(file, 'float', n=header$nBands)
        header$hiEdges <- pamBinRead(file, 'float', n=header$nBands)
    }
    return(header)
}