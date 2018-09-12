#' @title Read Clip Data
#' 
#' @description Reads binary data stored by the Clip Generator module.
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
#' @author Taiki Sakai \email{taiki.sakai@noaa.gov}
#' 
readClipData <- function(fid, fileInfo, data, debug=FALSE, ...) {
    error <- FALSE
    
    tryCatch({
        dataLength <- pamBinRead(fid, 'int32', n=1)
        if(dataLength==0) {
            return(list(data=data, error=error))
        }
        
        version <- fileInfo$moduleHeader$version
        
        if(version <= 1) {
            data$startSample <- pamBinRead(fid, 'int64', n=1)
            data$channelMap <- pamBinRead(fid, 'int32', n=1)
        }
        
        data$triggerMillis <- pamBinRead(fid, 'int64', n=1)
        
        if(version <= 1) {
            data$sampleDuration <- pamBinRead(fid, 'int32', n=1)
        }
        
        data$filename <- readJavaUTFString(fid)$str
        data$triggerName <- readJavaUTFString(fid)$str
        
        # Check if the object type = 2. If it is, there must be wav data at
        # the end of this object
        if(data$identifier == 2) {
            data$nChan <- pamBinRead(fid, 'int16', n=1)
            data$nSamps <- pamBinRead(fid, 'int32', n=1)
            data$scale <- 1 / pamBinRead(fid, 'float', n=1)
            data$wave <- matrix(
                pamBinRead(fid, 'int8', data$nSamps * data$nChan),
                nrow=data$nSamps, ncol=data$nChan) * data$scale
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