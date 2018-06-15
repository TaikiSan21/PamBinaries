#' @title Read Difar Data
#' 
#' @description Reads binary data stored by the Difar Processing module.
#'   
#' @param fid binary file identifier
#' @param fileInfo structure holding the file header and module header
#' @param data a structure containing standard data
#' @param skipLarge a flag of whether or not to skip reading the waveform
#' 
#' @return a structure containing data from a single object, and a logical
#'   flag if an error has occurred
#' 
#' @author Taiki Sakai \email{taiki.sakai@noaa.gov}
#' 
readDifarData <- function(fid, fileInfo, data, skipLarge = FALSE) {
    error <- FALSE
    
    tryCatch({
        dataLength <- pamBinRead(fid, 'int32', n=1)
        if(dataLength==0) {
            return(list(data=data, error=error))
        }
        
        version <- fileInfo$moduleHeader$version
        
        if(version <= 1) {
            data$startSample <- pamBinRead(fid, 'int64', n=1)
        }
        
        data$clipStart <- pamBinRead(fid, 'int64', n=1)
        
        if(version <= 1) {
            data$channelMap <- pamBinRead(fid, 'int32', n=1)
        }
        
        data$displaySampleRate <- pamBinRead(fid, 'float', n=1)
        data$demuxedLength <- pamBinRead(fid, 'int32', n=1)
        
        if(version <= 1) {
            minFreq <- pamBinRead(fid, 'float', n=1)
            maxFreq <- pamBinRead(fid, 'float', n=1)
            data$freqLimits <- c(minFreq, maxFreq)
        }
        
        data$amplitude <- pamBinRead(fid, 'float', n=1)
        data$gain <- pamBinRead(fid, 'float', n=1)
        data$selAngle <- pamBinRead(fid, 'float', n=1)
        data$selFreq <- pamBinRead(fid, 'float', n=1)
        data$speciesCode <- readJavaUTFString(fid)$str
        
        if(version >= 1) {
            data$trackedGroup <- readJavaUTFString(fid)$str
        }
        
        data$maxVal <- pamBinRead(fid, 'float', n=1)
        if(skipLarge) {
            seek(fid, data$demuxedLength * 3 * 2, origin = 'current')
        } else {
        data$demuxData <- matrix(
            pamBinRead(fid, 'int16', n=data$demuxedLength * 3),
            nrow=data$demuxedLength, ncol=3) * data$maxVal / 32767
        }
        
        data$numMatches <- pamBinRead(fid, 'int16', n=1)
        if(data$numMatches > 0) {
            data$latitude <- pamBinRead(fid, 'float', n=1)
            data$longitude <- pamBinRead(fid, 'float', n=1)
            
            if(fileInfo$moduleHeader$version >= 1) {
                errorX <- pamBinRead(fid, 'float', n=1)
                errorY <- pamBinRead(fid, 'float', n=1)
                data$errors <- c(errorX, errorY, 0)
            }
            
            for(i in 1:(data$numMatches-1)) { # This seems wrong, why -1?
                data$matchChan[i] <- pamBinRead(fid, 'int16', n=1)
                data$matchTime[i] <- pamBinRead(fid, 'int64', n=1)
            }
        }
        
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