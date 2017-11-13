#' @title Read Click Data
#' 
#' @description Reads binary data stored by the Click Detector module.
#'   
#' @param fid binary file identifier
#' @param fileInfo structure holding the file header and module header
#' @param data a structure containing standard data
#' @param getWave a flag for whether or not wave file data should be read
#' @param onlyWave a flag for skipping over other data when reading wave files
#' 
#' @return a structure containing data from a single object, and a logical
#'   flag if an error has occurred
#' 
#' @author Taiki Sakai \email{taiki.sakai@noaa.gov}
#' 
readClickData <- function(fid, fileInfo, data, getWave=FALSE, onlyWave=FALSE) {
    error <- FALSE
    tryCatch({
        dataLength <- pamBinRead(fid, 'int32', n=1, seek=onlyWave)
        if(dataLength==0) {
            return(list(data=data, error=error))
        }
        
        version <- fileInfo$moduleHeader$version
        
        if(version <= 3) {
            data$startSample <- pamBinRead(fid, 'int64', n=1, seek=onlyWave)
            data$channelMap <- pamBinRead(fid, 'int32', n=1)
        }
        
        data$triggerMap <- pamBinRead(fid, 'int32', n=1, seek=onlyWave)
        data$type <- pamBinRead(fid, 'int16', n=1, seek=onlyWave)
        
        if(version >= 2) {
            data$flags <- pamBinRead(fid, 'int32', n=1, seek=onlyWave)
        } else data$flags <- 0
        
        if(version <= 3) {
            nDelays <- pamBinRead(fid, 'int16', n=1)
            ################
            # matlab has if(nDelays), should be fine if nDelays is 0, so not needed
            ###################
            if(nDelays > 0) {
                data$delays <- pamBinRead(fid, 'float', n=nDelays, seek=onlyWave) #### THIS IS POSSIBLE ERROR NUMERIC - FLOAT
            }
        }
        
        nAngles <- pamBinRead(fid, 'int16', n=1)
        # if(nAngles) again
        if(nAngles > 0) {
            data$angles <- pamBinRead(fid, 'float', n=nAngles, seek=onlyWave)
        }
        
        if(version >= 3) {
            nAngleErrors <- pamBinRead(fid, 'int16', n=1)
            if(nAngleErrors > 0) {
                data$angleErrors <- pamBinRead(fid, 'float', n=nAngleErrors, seek=onlyWave)
            }
        } else data$angleErrors <- numeric() #unsure if equiv. to []
        
        if(version <= 3) {
            data$duration <- pamBinRead(fid, 'int16', n=1)
        } else data$duration <- data$sampleDuration
        
        data$nChan <- countChannels(data$channelMap)
        maxVal <- pamBinRead(fid, 'float', n=1)
        if(getWave) {
            data$wave <- matrix(
                pamBinRead(fid, 'int8', n = data$duration * data$nChan),
                nrow = data$duration, ncol = data$nChan) * maxVal / 127
        } else {
            seek(fid, data$duration * data$nChan, origin='current')
        }
        # Check if this matrix output is correct
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
