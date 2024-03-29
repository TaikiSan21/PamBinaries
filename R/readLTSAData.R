#' @title Read LTSA Data
#' 
#' @description Reads binary data stored by the LTSA module.
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
readLTSAData <- function(fid, fileInfo, data, debug=FALSE, ...) {
    error <- FALSE
    a <- 127*2/log(32767)
    b <- -127
    tryCatch({
        dataLength <- pamBinRead(fid, 'int32', n=1)
        if(dataLength==0) {
            return(list(data=data, error=error))
        }
        
        version <- fileInfo$moduleHeader$version
        
        if(version <= 1) {
            data$startSample <- pamBinRead(fid, 'int64', n=1)
        }
        
        if(version==0) {
            data$duration <- pamBinRead(fid, 'int64', n=1)
        }
        
        if(version <= 1) {
            data$channelMap <- pamBinRead(fid, 'int32', n=1)
        }
        
        data$endMillis <- pamBinRead(fid, 'int64', n=1)
        data$endDate <- millisToDateNum(data$endMillis)
        data$nFFT <- pamBinRead(fid, 'int32', n=1)
        data$maxVal <- pamBinRead(fid, 'float', n=1)
        
        # Version 0 scaled the data linearly to 16 bit
        if(version==0) {
            data$byteData <- pamBinRead(fid, 'int16', n = fileInfo$moduleHeader$fftLength / 2)
            data$data <- data$byteData / 32767 * data$maxVal
        } else {
            # After version 0, the data was first scaled to 16 bit and then
            # converted to a log so that it could be saved as an 8 bit
            data$byteData <- pamBinRead(fid, 'int8', n = fileInfo$moduleHeader$fftLength / 2)
            data$data <- exp((data$byteData - b)/a)* data$maxVal / 32767
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