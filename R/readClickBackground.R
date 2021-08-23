#' @title Read Click Detector Background data
#' 
#' @description Reads in the background data saved by the Click Detector
#'   
#' @param fid binary file identifier
#' @param fileInfo structure holding the file header, module header, and the
#'   appropriate function to read module specific data
#' @param data a structure containing standard data
#' 
#' @return a structure containing data from a single object
#' 
#' @author Michael Oswald \email{taiki.sakai@@noaa.gov}
#' 
readClickBackground <- function(fid, fileInfo, data) {
    
    # initialize variables
    error <- FALSE
    
    dataLength <- pamBinRead(fid, 'int32', n=1)
    data$noiseLen <- pamBinRead(fid, 'int16', n=1)
    data$backGround <- pamBinRead(fid, 'float', n=data$noiseLen)
    
    return(list(data=data, error=error))
}