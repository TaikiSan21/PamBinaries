#' @title Read Click Classifier Annotation
#' 
#' @description Reads binary data stored by Click Classifier annotations
#'   
#' @param fid binary file identifier
#' @param fileInfo structure holding the file header and module header
#' @param data a structure containing standard data
#' @param skipLarge a flag for whether or not to skip reading large wave file
#' @param getWave DEPRECATED: see skipLarge 
#' @param onlyWave DEPRECATED: see skipLarge
#' 
#' @return a structure containing data from a single object, and a logical
#'   flag if an error has occurred
#' 
#' @author Taiki Sakai \email{taiki.sakai@noaa.gov}
#' 
clickClassifierAnnotation <- function(fid, fileInfo, data) {
    error <- FALSE
    tryCatch({
        nClassifications <- pamBinRead(fid, 'int16', n=1)
        data$classifySet <- rep(0, nclassifications)
        for(i in 1:nClassifications) {
            data$classifySet[i] <- pamBinRead(fid, 'int16', n=1)
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