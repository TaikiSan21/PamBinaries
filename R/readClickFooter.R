#' @title Read Click Footer
#' 
#' @description Reads module footer information for the Click Detector module.
#'   Note that sometimes there is no additional footer information, so check
#'   first whether or not the binaryLength variable is 0.
#'   
#' @param file binary file to be read
#' 
#' @return footer information for Click Detector module
#' 
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#' 
readClickFooter <- function(file) {
    footer <- readStdModuleFooter(file)
    if(footer$binaryLength != 0) {
        footer$typesCountLength <- pamBinRead(file, 'int16', n=1)
        footer$typesCount <- pamBinRead(file, 'int32', n=footer$typesCountLength)
    }
    return(footer)
}