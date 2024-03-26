#' @title Add Frequency and Time to Pamguard Whistle Binaries
#' 
#' @description Adds items \code{freq} and \code{time} to a Pamguard binary file
#'   from the Whistle & Moan Detector
#'   
#' @param data either a \code{PamBinary} class object or just the \code{$data} from
#'   a PamBinary object
#' @param verbose logical flag to print calculated parameters
#' 
#' @return \code{data} with items \code{freq} and \code{time} added. These use the 
#'   calculated FFT window length, hope size, and sample rate to compute the frequency
#'   and time values of the saved whistle contour
#' 
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#' 
#' @examples 
#' 
#' # load example whistle file
#' wmFile <- system.file('extdata', 'WM.pgdf', package='PamBinaries')
#' wmData <- loadPamguardBinaryFile(wmFile)
#' # converts contour and FFT slice numbers to frequency and time values
#' wmData <- contourToFreq(wmData)
#' wmData$data[[1]]$contour
#' wmData$data[[1]]$freq
#' wmData$data[[1]]$time
#' 
#' @importFrom dplyr bind_rows distinct
#' @export
#' 
contourToFreq <- function(data, verbose=FALSE) {
    if(inherits(data, 'PamBinary')) {
        data$data <- contourToFreq(data$data, verbose)
        return(data)
    }
    if(length(data) == 0) {
        return(data)
    }
    if(!(all(c('sliceData', 'nSlices', 'sampleDuration', 'startSample', 'maxFreq') %in%
             names(data[[1]])))) {
        stop('Appears data is not a Whistle and Moan Detector binary file.')
    }
    fftParams <- getPamFft(data, method='new')
    if(is.null(fftParams)) {
        fftParams <- list(sr=NA, fftLen=NA, fftHop=NA)
    }
    sr <- fftParams$sr
    fftLen <- fftParams$wl
    fftHop <- fftParams$hop
    if(verbose) {
        cat('SR: ', sr, ' Len: ', fftLen, ' Hop: ', fftHop, sep ='')
    }
    for(i in seq_along(data)) {
        data[[i]]$freq <- data[[i]]$contour * sr / fftLen
        data[[i]]$allFreq <- do.call(cbind, lapply(data[[i]]$sliceData, function(x) x$peakData)) * sr / fftLen
        data[[i]]$time <- sapply(data[[i]]$sliceData,
                                 function(x) x$sliceNumber) * fftHop / sr
    }
    data
}

getPamFft <- function(data, method=c('new', 'old')) {
    if(inherits(data, 'PamBinary')) {
        # data$data <- contourToFreq(data$data)
        return(getPamFft(data$data, method=method))
    }
    if(length(data) == 0) {
        return(NULL)
    }
    if(!(all(c('sliceData', 'nSlices', 'sampleDuration', 'startSample', 'maxFreq') %in%
             names(data[[1]])))) {
        # stop('Appears data is not a Whistle and Moan Detector binary file.')
        return(NULL)
    }
    if(length(data) == 1) {
        method <- 'old'
    }
    switch(
        match.arg(method),
        'new' = {
            samplePairs <- bind_rows(
                lapply(data, function(x) {
                    list(start=x$startSample+1, slice=x$sliceData[[1]]$sliceNumber)
                }))
            samplePairs <- distinct(samplePairs)
            if(nrow(samplePairs) <= 1) {
                return(getPamFft(data, method='old'))
            }
            
            hops <- diff(samplePairs$start) / diff(samplePairs$slice)
            hops <- unique(hops)
            hops <- hops[is.finite(hops)]
            if(length(hops) > 1) {
                hops <- as.integer(hops)
                checkPow2 <- round(log2(hops)) == log2(hops)
                if(!any(checkPow2)) {
                    fftHop <- hops[1]
                } else {
                    fftHop <- which(checkPow2)[1]
                }
            } else {
                fftHop <- hops
            }
            fftLen <- data[[1]]$sampleDuration -
                (data[[1]]$sliceData[[data[[1]]$nSlices]]$sliceNumber - data[[1]]$sliceData[[1]]$sliceNumber) * fftHop
            sr <- fftLen * data[[1]]$maxFreq /
                max(unlist(lapply(data[[1]]$sliceData, function(x) x$peakData)))
            if(any(c(fftHop, fftLen, sr) <= 0) ||
               any(is.na(c(fftHop, fftLen, sr))) ||
               any(is.infinite(c(fftHop, fftLen, sr)))) {
                return(getPamFft(data, method='old'))
            }
            # return(list(sr=sr, hop=fftHop, wl=fftLen))
            
        },
        'old' = {
            tempData <- data[[1]]
            if(tempData$sliceData[[1]]$sliceNumber == 0) {
                if(length(data) == 1) {
                    return(NULL)
                }
                tempData <- data[[2]]
            }
            fftHop <- (tempData$startSample + 1)/tempData$sliceData[[1]]$sliceNumber
            fftLen <- tempData$sampleDuration -
                (tempData$sliceData[[tempData$nSlices]]$sliceNumber - tempData$sliceData[[1]]$sliceNumber) * fftHop
            sr <- fftLen * tempData$maxFreq /
                max(unlist(lapply(tempData$sliceData, function(x) x$peakData)))
            # return(list(sr=sr, hop=fftHop, wl=fftLen))
            if(any(c(fftHop, fftLen, sr) <= 0) ||
               any(is.na(c(fftHop, fftLen, sr))) ||
               any(is.infinite(c(fftHop, fftLen, sr)))) {
                return(NULL)
            }
        }
    )
    list(sr=sr, hop=fftHop, wl=fftLen)
}
