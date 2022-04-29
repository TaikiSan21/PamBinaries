#' @title Read Beam Former Annotation
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
readBeamFormerAnnotation <- function(fid, fileInfo, anVersion, debug=FALSE, ...) {
    error <- FALSE
    result <- list()
    tryCatch({
        result$hydrophones <- pamBinRead(fid, 'uint32', n=1)
        result$arrayType <- pamBinRead(fid, 'int16', n=1)
        result$localisationContent <- pamBinRead(fid, 'uint32', n=1)
        result$nAngles <- pamBinRead(fid, 'int16', n=1)
        result$angles <- pamBinRead(fid, 'float', n=result$nAngles)
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
