#' @title Read TDBL Annotation
#' 
#' @description Reads binary data stored by TDBL annotation module
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
readTDBLAnnotation <- function(fid, fileInfo, anVersion, debug=FALSE, ...) {
    error <- FALSE
    result <- list()
    tryCatch({
        nAngles <- pamBinRead(fid, 'int16', n=1)
        result$angles <- pamBinRead(fid, 'float', n=nAngles)
        nErrors <- pamBinRead(fid, 'int16', n=1)
        result$angleErrors <- pamBinRead(fid, 'float', n=nErrors)
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
