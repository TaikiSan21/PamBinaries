#' @title Read Noise Monitor Data
#' 
#' @description Reads binary data stored by the Noise Monitor.
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
readNoiseMonData <- function(fid, fileInfo, data) {
    error <- FALSE
    
    tryCatch({
        dataLength <- pamBinRead(fid, 'int32', n=1)
        if(dataLength==0) {
            return(list(data=data, error=error))
        }
        
        version <- fileInfo$moduleHeader$version
        
        data$iChan <- pamBinRead(fid, 'int16', n=1)
        data$nBands <- pamBinRead(fid, 'int16', n=1)
        
        if(version >= 1) {
            data$nMeasures <- pamBinRead(fid, 'int16', n=1)
        } else {
            data$nMeasures <- 4
        }
        
        if(version <= 1) {
            n <- pamBinRead(fid, 'float', n = data$nBands * data$nMeasures)
        } else {
            n <- pamBinRead(fid, 'int16', n = data$nBands * data$nMeasures) / 100
        }
        
        data$noise <- matrix(n, nrow = data$nMeasures, ncol = data$nBands)
        
        return(list(data=data, error=error))
    # }, warning = function(w) {
    #     print(paste('Warning occurred: ', w))
    #     return(list(data=data, error=error))
    }, error = function(e) {
        print(paste('Error reading ', fileInfo$fileHeader$moduleType, ' data object. Data read:'))
        print(data)
        print(e)
        error <- TRUE
        return(list(data=data, error=error))
    })
}