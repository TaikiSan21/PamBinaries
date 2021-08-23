#' @title Read Spectral Background data
#' 
#' @description Reads in the background data saved by various detectors (e.g
#'    WMD, Right Whale Edge Detector, etc) EXCEPT FOR the Click Detector
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
readSpectralBackground <- function(fid, fileInfo, data) {
    
    # initialize variables
    error <- FALSE
    
    dataLength <- pamBinRead(fid, 'int32', n=1)
    data$firstBin <- pamBinRead(fid, 'int32', n=1)
    data$noiseLen <- pamBinRead(fid, 'int32', n=1)
    data$backGround <- pamBinRead(fid, 'float', n=data$noiseLen)
    
    return(list(data=data, error=error))
}