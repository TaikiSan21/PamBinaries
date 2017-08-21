#' @title Read AIS Data
#' 
#' @description Reads binary data stored by the AIS Processing module.
#'   
#' @param fid binary file identifier
#' @param fileInfo structure holding the file header and module header
#' @param data a structure containing standard data
#' 
#' @return a structure containing data from a single object, and a logical
#'   flag if an error has occurred
#' 
#' @author Taiki Sakai \email{taiki.sakai@noaa.gov}
#' 
readAISData <- function(fid, fileInfo, data) {
    error <- FALSE
    
    tryCatch({
        # Read AIS Processing Module specific data
        dataLength <- pamBinRead(fid, 'int32', n=1)
        
        if(dataLength==0) {
            return(list(data=data, error=error))
        }
        
        data$mmsiNumber <- pamBinRead(fid, 'int32', n=1)
        data$fillBits <-pamBinRead(fid, 'int16', n=1)
        data$charData <- readJavaUTFString(fid)$str
        data$aisChannel <- readJavaUTFString(fid)$str
        return(list(data=data, error=error))
    }, warning = function(w) {
        print(paste('Warning occurred: ', w))
        return(list(data=data, error=error))
    }, error = function(e) {
        print('Error reading AIS data object. Data read:')
        print(data)
        print(e)
        error <- TRUE
        return(list(data=data, error=error))
    })
}
