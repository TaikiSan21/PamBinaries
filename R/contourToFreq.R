#' @title Add Frequency and Time to Pamguard Whistle Binaries
#' 
#' @description Adds items \code{freq} and \code{time} to a Pamguard binary file
#'   from the Whistle & Moan Detector
#'   
#' @param data either a \code{PamBinary} class object or just the \code{$data} from
#'   a PamBinary object
#' 
#' @return \code{data} with items \code{freq} and \code{time} added. These use the 
#'   calculated FFT window length, hope size, and sample rate to compute the frequency
#'   and time values of the saved whistle contour
#' 
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#' 
#' @export
#' 
contourToFreq <- function(data) {
    if(inherits(data, 'PamBinary')) {
        data$data <- contourToFreq(data$data)
        return(data)
    }
    if(length(data) == 0) {
        return(data)
    }
    if(!(all(c('sliceData', 'nSlices', 'sampleDuration', 'startSample', 'maxFreq') %in%
             names(data[[1]])))) {
        stop('Appears data is not a Whistle and Moan Detector binary file.')
    }
    tempData <- data[[1]]
    if(tempData$sliceData[[1]]$sliceNumber == 0) {
        tempData <- data[[2]]
    }
    fftHop <- (tempData$startSample + 1)/tempData$sliceData[[1]]$sliceNumber
    fftLen <- tempData$sampleDuration - (tempData$nSlices - 1) * fftHop
    sr <- fftLen * tempData$maxFreq /
        max(unlist(lapply(tempData$sliceData, function(x) x$peakData)))
    for(i in seq_along(data)) {
        data[[i]]$freq <- data[[i]]$contour * sr / fftLen
        data[[i]]$allFreq <- do.call(cbind, lapply(data[[i]]$sliceData, function(x) x$peakData)) * sr / fftLen
        data[[i]]$time <- sapply(data[[i]]$sliceData,
                                 function(x) x$sliceNumber) * fftHop / sr
    }
    data
}
