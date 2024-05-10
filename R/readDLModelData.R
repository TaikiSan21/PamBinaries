#' @title Read Deep Learning Model Data
#' 
#' @description Reads binary data stored by the Deep Learning Model module
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
readDLModelData <- function(fid, fileInfo, data, debug=FALSE, ...) {
    error <- FALSE
    
    tryCatch({
        dataLength <- pamBinRead(fid, 'int32', n=1)
        if(dataLength==0) {
            return(list(data=data, error=error))
        }
        data$tpye <- pamBinRead(fid, 'int8', n=1)
        data$isbinary <- pamBinRead(fid, 'int8', n=1) > 0
        scale <- pamBinRead(fid, 'float32', n=1)
        nSpecies <- pamBinRead(fid, 'int16', n=1)
        pred <- vector('numeric', length=nSpecies)
        for(i in 1:nSpecies) {
            pred[i] <- pamBinRead(fid, 'int16', n=1)/scales
        }
        data$predictions <- pred
        nclass <- pamBinRead(fid, 'int16', n=1)
        if(nclass > 0) {
            for(i in 1:nclass) {
                pamBinRead(fid, 'int16', n=1)
            }
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
