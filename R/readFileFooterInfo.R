#' @title Read File Footer
#' 
#' @description Reads in the binary file footer. The input variable version
#'   is the file format read in from the file header. As of version 3, the 
#'   file footer includes the lowest and highest UID values in the file.
#'   
#' @param fid binary file to be read
#' @param version binary file version
#' 
#' @return footer information common to all files
#' 
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#' 
readFileFooterInfo <- function(fid, version) {
    footer <- list()
    footer$length <- pamBinRead(fid, 'int32', n=1)
    footer$identifier <- pamBinRead(fid, 'int32', n=1)
    footer$nObjects <- pamBinRead(fid, 'int32', n=1)
    footer$dataDate <- millisToDateNum(pamBinRead(fid, 'int64', n=1))
    footer$analysisDate <- millisToDateNum(pamBinRead(fid, 'int64', n=1))
    footer$endSample <- pamBinRead(fid, 'int64', n=1)
    if(version >= 3) {
        footer$lowestUID <- pamBinRead(fid, 'int64', n=1)
        footer$highestUID <- pamBinRead(fid, 'int64', n=1)
    }
    footer$fileLength <- pamBinRead(fid, 'int64', n=1)
    footer$endReason <- pamBinRead(fid, 'int32', n=1)
    return(footer)
}