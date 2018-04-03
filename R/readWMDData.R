#' @title Read Whistle and Moan Data
#' 
#' @description Reads binary data stored by the Whistle & Moan Detector
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
readWMDData <- function(fid, fileInfo, data) {
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
        
        data$nSlices <- pamBinRead(fid, 'int16', n=1)
        
        if(version >= 1) {
            data$amplitude <- pamBinRead(fid, 'int16', n=1) / 100
        }
        
        if(version == 1) {
            data$nDelays <- pamBinRead(fid, 'int8', n=1)
            data$delays <- pamBinRead(fid, 'int16', n = data$nDelays) # need to scale this still !!!
            data$delays <- data$delays / fileInfo$moduleHeader$delayScale # Not sure if this is right - diff from matlab
        }
        
        data$sliceData <- list()
        data$contour <- rep(0, data$nSlices)
        data$contWidth <- rep(0, data$nSlices)
        for(i in 1:data$nSlices) {
            aSlice <- list()
            aSlice$sliceNumber <- pamBinRead(fid, 'int32', n=1)
            aSlice$nPeaks <- pamBinRead(fid, 'int8', n=1)
            aSlice$peakData <- matrix(0, nrow=4, ncol=aSlice$nPeaks)
            for(p in 1:aSlice$nPeaks) {
                sss <- pamBinRead(fid, 'int16', n=4)
                aSlice$peakData[,p] <- sss
            }
            data$sliceData[[i]] <- aSlice
            data$contour[i] <- aSlice$peakData[2,1]
            data$contWidth[i] <- aSlice$peakData[3,1] - aSlice$peakData[1,1] + 1
        }
        data$meanWidth <- mean(data$contWidth)
        
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