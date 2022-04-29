#' @title Read Deep Learning Detection Data
#' 
#' @description Reads binary data stored by the DbHt module.
#'   
#' @param fid binary file identifier
#' @param fileInfo structure holding the file header and module header
#' @param data a structure containing standard data
#' @param debug logical flag to show more info on errors
#' @param \dots Arguments passed to other functions
#' 
#' @return a structure containing data from a single object, and a logical
#'   flag if an error has occurred
#' 
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#' 
readDLDetData <- function(fid, fileInfo, data, debug=FALSE, ...) {
    error <- FALSE
    
    tryCatch({
        dataLength <- pamBinRead(fid, 'int32', n=1)
        if(dataLength==0) {
            return(list(data=data, error=error))
        }
        data$nChan <- pamBinRead(fid, 'int16', n=1)
        data$nSamps <- pamBinRead(fid, 'int32', n=1)
        data$scale <- 1/pamBinRead(fid, 'float', n=1)
        data$wave <- matrix(pamBinRead(fid, 'int8', n=data$nSamps*data$nChan), ncol=data$nChan) / data$scale
        
        return(list(data=data, error=error))
        # }, warning = function(w) {
        #     print(paste('Warning occurred: ', w))
        #     return(list(data=data, error=error))
    }, error = function(e) {
        if(debug) {
            print(paste0('Error reading ', fileInfo$fileHeader$moduleType, ' Data read:'))
            print(data)
            print(e)
        }
        error <- TRUE
        return(list(data=data, error=error))
    })
}
