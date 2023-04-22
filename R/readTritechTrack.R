#' @title Read a Tritech Track
#' 
#' @description Reads binary data stored by the Gemini Tritech Module
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
#' @author Benjamin Blundell \email{bjb8@st-andrews.ac.uk}
#' 
readTritechTrack <- function(fid, fileInfo, data, debug=FALSE, ...) {
    error <- FALSE
    
    tryCatch({
        dataLength <- pamBinRead(fid, 'int32', n=1)
        if(dataLength==0) {
            return(list(data=data, error=error))
        }
        
        version <- fileInfo$moduleHeader$version

        data$nPoints <- pamBinRead(fid, 'int32', n=1)
        data$nSonar <- pamBinRead(fid, 'int8', n=1)
        data$sonarIds <- pamBinRead(fid, 'int16', n=data$nSonar)
        data$straightLength <- pamBinRead(fid, 'float', n=1)
        data$wobblyLength <- pamBinRead(fid, 'float', n=1)
        data$meanOccupancy <- pamBinRead(fid, 'float', n=1)

        data$track <- data.frame(timeMillis = rep(0, data$nPoints),
                                sonarId = rep(0, data$nPoints),
                                minBearing = rep(0, data$nPoints),
                                maxBearing = rep(0, data$nPoints),
                                peakBearing = rep(0, data$nPoints),
                                minRange = rep(0, data$nPoints),
                                maxRange = rep(0, data$nPoints),
                                peakRange = rep(0, data$nPoints),
                                objSize = rep(0, data$nPoints),
                                occupancy = rep(0, data$nPoints),
                                aveValue = rep(0, data$nPoints),
                                totValue = rep(0, data$nPoints),
                                maxValue = rep(0, data$nPoints))

        for(i in 1:data$nPoints) {
            data$track$timeMillis[i] <- pamBinRead(fid, 'double', n=1)
            data$track$sonarId[i] <- pamBinRead(fid, 'int16', n=1)
            data$track$minBearing[i] <- pamBinRead(fid, 'float', n=1)
            data$track$maxBearing[i] <- pamBinRead(fid, 'float', n=1)
            data$track$peakBearing[i] <- pamBinRead(fid, 'float', n=1)
            data$track$minRange[i] <- pamBinRead(fid, 'float', n=1)
            data$track$maxRange[i] <- pamBinRead(fid, 'float', n=1)
            data$track$peakRange[i] <- pamBinRead(fid, 'float', n=1)
            data$track$objSize[i] <- pamBinRead(fid, 'float', n=1)
            data$track$occupancy[i] <- pamBinRead(fid, 'float', n=1)
            data$track$aveValue[i] <- pamBinRead(fid, 'int16', n=1)
            data$track$totValue[i] <- pamBinRead(fid, 'int32', n=1)
            data$track$maxValue[i] <- pamBinRead(fid, 'int16', n=1)
        }

        data$date <- millisToDateNum(data$millis)
         
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