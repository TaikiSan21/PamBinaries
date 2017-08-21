#' @title Read Standard Module Header
#' 
#' @description Reads the module header information common to all modules.
#'   Differs from the legacy code in that it does not read in or skip any
#'   information specific to a module.
#'   
#' @param file binary file to be read
#' 
#' @return header information common to all modules
#' 
#' @author Taiki Sakai \email{taiki.sakai@noaa.gov}
#' 
readStdModuleHeader <- function(file) {
    header <- list()
    header$length <- pamBinRead(file, 'int32', n=1)
    header$identifier <- pamBinRead(file, 'int32', n=1)
    header$version <- pamBinRead(file, 'int32', n=1)
    header$binaryLength <- pamBinRead(file, 'int32', n=1)
    return(header)
}