#' @title Read Deep Learning Annotation
#' 
#' @description Reads binary data stored by the DbHt module.
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
readDLAnnotation <- function(fid, fileInfo, anVersion, debug=FALSE, ...) {
    error <- FALSE
    result <- list()
    tryCatch({
        numModels <- pamBinRead(fid, 'int16', n=1)
        for(i in 1:numModels) {
            result[[i]] <- readModelData(fid)
        }
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

readModelData <- function(fid) {
    result <- list()
    # original matlab has 'char' here but that doesnt seem right?
    # all comparisons after are to int 0/1/2
    # modelType <- pamBinRead(fid, 'char', n=1)
    # isBinary <- pamBinRead(fid, 'char', n=1)
    modelType <- pamBinRead(fid, 'int8', n=1)
    isBinary <- pamBinRead(fid, 'int8', n=1)
    isBinary <- isBinary != 0
    scale <- pamBinRead(fid, 'float', n=1)
    nSpecies <- pamBinRead(fid, 'int16', n=1)
    
    data <- rep(-1, nSpecies)
    for(i in 1:nSpecies) {
        data[i] <- pamBinRead(fid, 'int16', n=1) / scale
    }
    
    nClass <- pamBinRead(fid, 'int16', n=1)
    classNames <- rep(-1, nClass)
    for(i in 1:nClass) {
        classNames[i] <- pamBinRead(fid, 'int16', n=1)
    }
    
    switch(as.character(modelType),
           '0' = { # generic ddep learning anno
               result$predictions <- data
               result$classID <- classNames
               result$isBinary <- isBinary
               result$type <- modelType
           },
           '1' = { #sound spot classifier
               result$predictions <- data
               result$classID <- classNames
               result$isBinary <- isBinary
               result$type <- modelType
           },
           '2' = { #dummy result
               result$predictions <- numeric(0)
               result$type <- 'dummy'
           }
    )
    result
}