#' @title Read GPL Detections
#' 
#' @description Reads binary data stored by the GPL Module.
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
#' @author Michael Oswald \email{mo55@@st-andrews.ac.uk}
#' 
readGPLDetections <- function(fid, fileInfo, data, debug=FALSE, ...) {
    error <- FALSE
    
    tryCatch({
        skipDummy <- pamBinRead(fid, 'int32', n=1)
        dataLength <- pamBinRead(fid, 'int32', n=1)
        if(dataLength==0) {
            return(list(data=data, error=error))
        }
        
        data$timeRes <- pamBinRead(fid, 'float', n=1)
        data$freqRes <- pamBinRead(fid, 'float', n=1)
        data$area <- pamBinRead(fid, 'int16', n=1)
        bitDepth <- pamBinRead(fid, 'int8', n=1)
        
        if(bitDepth == 8) {
            pType <- 'uint8'
        }else {
            pType <- 'uint16'
        }
        
        data$energy <- rep(0, data$area)
        data$points <- matrix(0, nrow=2, ncol=data$area)
        
        for(i in 1:data$area) {
            sss <- pamBinRead(fid, pType, n=2)
            data$points[,i] <- sss
            data$energy[i] <- pamBinRead(fid, 'float', n=1)
        }
        
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
