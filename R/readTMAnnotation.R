#' @title Read Target Motion Annotation
#' 
#' @description Reads binary data stored by beam former annotation module
#'   
#' @param fid binary file identifier
#' @param fileInfo structure holding the file header and module header
#' @param anVersion annotation version
#' @param debug logical flag to show more info on errors
#' @param \dots Arguments passed to other functions
#' 
#' @return a structure containing data from a single object, and a logical
#'   flag if an error has occurred
#' 
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#' 
readTMAnnotation <- function(fid, fileInfo, anVersion, debug=FALSE, ...) {
    error <- FALSE
    result <- list()
    tryCatch({
        result$model <- readJavaUTFString(fid)
        result$nLocations <- pamBinRead(fid, 'int16', n=1)
        result$hydrophones <- pamBinRead(fid, 'uint32', n=1)
        loc <- list(latitude=numeric(0),
                    longitude=numeric(0),
                    height=numeric(0),
                    error= character(0))
        for(i in 1:result$nLocations) {
            loc$latitude[i] <- pamBinRead(fid, 'double', n=1)
            loc$longitude[i] <- pamBinRead(fid, 'double', n=1)
            loc$height[i] <- pamBinRead(fid, 'float', n=1)
            loc$error[i] <- readJavaUTFString(fid)
        }
        result$loc <- loc
        return(result)
    }, error = function(e) {
        if(debug) {
            print(paste0('Error reading ', fileInfo$fileHeader$moduleType, ' Data read:'))
            print(result)
            print(e)
        }
        error <- TRUE
        return(result)
    })
}
