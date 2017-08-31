#' @title Read Click Data
#' 
#' @description Reads binary data stored by the Click Detector module.
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
readClickData <- function(fid, fileInfo, data) {
    error <- FALSE
    tryCatch({
        dataLength <- pamBinRead(fid, 'int32', n=1)
        if(dataLength==0) {
            return(list(data=data, error=error))
        }
        
        version <- fileInfo$moduleHeader$version
        
        if(version <= 3) {
            data$startSample <- pamBinRead(fid, 'int64', n=1)
            data$channelMap <- pamBinRead(fid, 'int32', n=1)
        }
        
        data$triggerMap <- pamBinRead(fid, 'int32', n=1)
        data$type <- pamBinRead(fid, 'int16', n=1)
        
        if(version >= 2) {
            data$flags <- pamBinRead(fid, 'int32', n=1)
        } else data$flags <- 0
        
        if(version <= 3) {
            nDelays <- pamBinRead(fid, 'int16', n=1)
            ################
            # matlab has if(nDelays), should be fine if nDelays is 0, so not needed
            ###################
            data$delays <- pamBinRead(fid, 'float', n=nDelays) #### THIS IS POSSIBLE ERROR NUMERIC - FLOAT
        }
        
        nAngles <- pamBinRead(fid, 'int16', n=1)
        # if(nAngles) again
        data$angles <- pamBinRead(fid, 'float', n=nAngles)
        
        if(version >= 3) {
            nAngleErrors <- pamBinRead(fid, 'int16', n=1)
            data$angleErrors <- pamBinRead(fid, 'float', n=nAngleErrors)
        } else data$angleErrors <- numeric() #unsure if equiv. to []
        
        if(version <= 3) {
            data$duration <- pamBinRead(fid, 'int16', n=1)
        } else data$duration <- data$sampleDuration
        
        data$nChan <- countChannels(data$channelMap)
        maxVal <- pamBinRead(fid, 'float', n=1)
        data$wave <- matrix(
            pamBinRead(fid, 'int8', n = data$duration * data$nChan),
            nrow = data$duration, ncol = data$nChan) * maxVal / 127 # Check if this matrix output is correct
        return(list(data=data, error=error))
    # }, warning = function(w) {
        # print(paste('Warning occurred: ', w))
        # return(list(data=data, error=error))
    }, error = function(e) {
        print('Error reading ', fileInfo$fileHeader$moduleType, ' Data read:')
        print(data)
        print(e)
        error <- TRUE
        return(list(data=data, error=error))
    })
}
