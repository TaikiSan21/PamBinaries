#' @title Read Matched Classifier Annotation
#' 
#' @description Reads annotations from the matched click classifier. The matched
#'    matched click classifier annotates click detections with a threshold, matchcorr
#'    and rejectcorr values. The threshold value is used in the binary classification
#'    process. If it exceeds a hard value then the click is classified with the
#'    set type. The matchcorr and rejectcorr values are simply the correlation 
#'    values of the match and reject templates with the click.
#'   
#' @param fid binary file identifier
#' @param fileInfo structure holding the file header and module header
#' @param anVersion version id of annotation module
#' @param debug logical flag to show more info on errors
#' 
#' @return a vector with the threshold, matchcorr, and rejectcorr values. See description.
#' 
#' @author Taiki Sakai \email{taiki.sakai@noaa.gov}
#' 
readMatchClsfrAnnotation <- function(fid, fileInfo, anVersion, debug=FALSE) {
    error <- FALSE
    data <- c()
    tryCatch({
        if(anVersion == 1) {
            threshold <- pamBinRead(fid, 'double', n=1)
            matchcorr <- pamBinRead(fid, 'double', n=1)
            rejectcorr <- pamBinRead(fid, 'double', n=1)
            data <- list('threshold'=threshold, 'matchcorr'=matchcorr, 'rejectcorr'=rejectcorr)
            return(data)
        }
        if(anVersion == 2) {
            nTemplates <- pamBinRead(fid, 'int16', n=1)
            data <- data.frame(threshold = rep(0, nTemplates),
                               mathccorr = rep(0, nTemplates),
                               rejectcorr = rep(0, nTemplates))
            for(i in 1:nTemplates) {
                data$threshold[i] <- pamBinRead(fid, 'double', n=1)
                data$matchcorr[i] <- pamBinRead(fid, 'double', n=1)
                data$rejectcorr[i] <- pamBinRead(fid, 'double', n=1)
            }
            return(data)
        }
    }, error = function(e) {
        if(debug) {
            print(paste0('Error reading ', fileInfo$fileHeader$moduleType, 
                         ' matched classifier annotation. Data read:'))
            print(data)
            print(e)
        }
        error <- TRUE
        return(data)
        # return(list(data=data, error=error))
    })
}