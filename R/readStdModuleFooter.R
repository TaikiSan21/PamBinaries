#' @title Read Standard Module Footer
#' 
#' @description Reads the module footer information common to all modules.
#'   Differs from the legacy code in that it does not read in or skip any
#'   information specific to a module.
#'   
#' @param file binary file to be read
#' 
#' @return footer information common to all modules
#' 
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#' 
readStdModuleFooter <- function(file) {
    footer <- list()
    footer$length <- pamBinRead(file, 'int32', n=1)
    footer$identifier <- pamBinRead(file, 'int32', n=1)
    footer$binaryLength <- pamBinRead(file, 'int32', n=1)
    return(footer)
}