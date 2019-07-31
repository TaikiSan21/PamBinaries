#' @title Read Click Trigger Level
#' 
#' @description Reads binary data stored by the click detector trigger
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
readClickTriggerData <- function(fid, fileInfo, data, debug=FALSE, ...) {
    error <- FALSE
    nChan <- countChannels(data$channelMap)
    tryCatch({
        dataLength <- pamBinRead(fid, 'int32', n=1)
        if(dataLength==0) {
            return(list(data=data, error=error))
        }
        scale <- pamBinRead(fid, 'float', n=1)
        data$rawLevels <- pamBinRead(fid, 'int16', n=nChan) / scale
        cal <- fileInfo$moduleHeader$calibration
        if(!is.null(cal)) {
            data$absLevelsdB <- 20*log10(data$rawLevels) + cal
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