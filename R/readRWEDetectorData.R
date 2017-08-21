#' @title Read Right Whale Edge Detector Data
#' 
#' @description Reads binary data stored by the Right Whale Edge Detector.
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
readRWEDetectorData <- function(fid, fileInfo, data) {
    error <- FALSE
    
    tryCatch({
        dataLength <- pamBinRead(fid, 'int32', n=1)
        if(dataLength==0) {
            return(list(data=data, error=error))
        }
        
        version <- fileInfo$moduleHeader$version
        
        if(version==0) {
            data$startSample <- pamBinRead(fid, 'int64', n=1)
            data$channelMap <- pamBinRead(fid, 'int32', n=1)
        }
        
        data$type <- pamBinRead(fid, 'int16', n=1)
        data$signal <- pamBinRead(fid, 'float', n=1)
        data$noise <- pamBinRead(fid, 'float', n=1)
        data$nSlices <- pamBinRead(fid, 'int16', n=1)
        
        zeros <- rep(0, data$nSlices)
        
        data$times <- zeros
        data$sliceNums <- zeros
        data$loFreqs <- zeros
        data$peakFreqs <- zeros
        data$hiFreqs <- zeros
        data$peakAmps <- zeros
        
        for(i in 1:data$nSlices) {
            data$sliceNums[i] <- pamBinRead(fid, 'int16', n=1)
            data$loFreqs[i] <- pamBinRead(fid, 'int16', n=1)
            data$peakFreqs[i] <- pamBinRead(fid, 'int16', n=1)
            data$hiFreqs[i] <- pamBinRead(fid, 'int16', n=1)
            data$peakAmps[i] <- pamBinRead(fid, 'float', n=1)
        }
        
        return(list(data=data, error=error))
    }, warning = function(w) {
        print(paste('Warning occurred: ', w))
        return(list(data=data, error=error))
    }, error = function(e) {
        print(paste('Error reading ', fileInfo$fileHeader$moduleType, ' data object. Data read:'))
        print(data)
        print(e)
        error <- TRUE
        return(list(data=data, error=error))
    })
}