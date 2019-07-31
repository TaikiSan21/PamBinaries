#' @title Read Click Classifier Annotation
#' 
#' @description Reads binary data stored by Click Classifier annotations
#'   
#' @param fid binary file identifier
#' @param fileInfo structure holding the file header and module header
#' @param debug logical flag to show more info on errors
#' 
#' @return a vector of click classifiers, represented by the click type flag
#' 
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#' 
readClickClsfrAnnotation <- function(fid, fileInfo, debug=FALSE) {
    error <- FALSE
    tryCatch({
        nClassifications <- pamBinRead(fid, 'int16', n=1)
        classifySet <- rep(0, nClassifications)
        for(i in 1:nClassifications) {
            classifySet[i] <- pamBinRead(fid, 'int16', n=1)
        }
        return(classifySet)
        # return(list(data=data, error=error))
        # }, warning = function(w) {
        #     print(paste('Warning occurred: ', w))
        #     return(list(data=data, error=error))
    }, error = function(e) {
        if(debug) {
            print(paste0('Error reading ', fileInfo$fileHeader$moduleType, ' classification set. Data read:'))
            print(classifySet)
            print(e)
        }
        error <- TRUE
        return(classifySet)
        # return(list(data=data, error=error))
    })
}