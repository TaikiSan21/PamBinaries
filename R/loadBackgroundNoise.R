#' @title Load and Format Background Noise Data
#' 
#' @description Reads and formats background noise data from Pamguard binary files
#'   or if not present in the original file will try to read the accompanying 
#'   .pgnf noise file if it exists
#' 
#' @param x character pointing to a Pamguard binary file, or a \code{PamBinary}
#'   object created by \link{loadPamguardBinaryFile}. For plotting, either of these
#'   or the output from \link{loadBackgroundNoise}
#'   
#' @return A list with \code{times} storing the POSIXct time of each background 
#'   measurement, and \code{background} a matrix of background values. For binary
#'   data based on spectrogram measurements, there will also be \code{freq} the 
#'   frequency in Hertz for each column of background measurement
#' 
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#' 
#' @examples 
#' 
#' # load the example click binary data, leaving date as numeric
#' gplFile <- system.file('extdata', 'GPL.pgdf', package='PamBinaries')
#' gplNoise <- loadBackgroundNoise(gplFile)
#' plotBackgroundNoise(gplNoise)
#' 
#' @importFrom graphics axis image lines
#' @export
#'
loadBackgroundNoise <- function(x) {
    if(is.character(x) &&
       file.exists(x)) {
        x <- loadPamguardBinaryFile(x)
    }
    info <- x$fileInfo
    oneDat <- x$data[[1]]
    if(is.null(info$background)) {
        noiseFile <- gsub('df$', 'nf', info$fileName)
        if(file.exists(noiseFile)) {
            noise <- loadPamguardBinaryFile(noiseFile)
            noise$data <- list(oneDat)
            return(loadBackgroundNoise(noise))
        }
        return(NULL)
    }
    type <- info$fileHeader$moduleType
    times <- convertPgDate(sapply(info$background, function(b) {
        b$date
    }))
    bgData <- matrix(NA, nrow = length(info$background), ncol = info$background[[1]]$noiseLen)
    for(i in seq_along(info$background)) {
        bgData[i, ] <- info$background[[i]]$backGround
    }
    result <- list(times = times, background=bgData)
    if(type == 'Click Detector') {
        return(result)
    }
    switch(type,
           'GPL Detector' = {
               freqRes <- oneDat$freqRes
           },
           'WhistlesMoans' = {
               if(oneDat$sliceData[[1]]$sliceNumber == 0) {
                   oneDat <- x$data[[2]]
               }
               fftHop <- (oneDat$startSample + 1)/oneDat$sliceData[[1]]$sliceNumber
               fftLen <- oneDat$sampleDuration -
                   (oneDat$sliceData[[oneDat$nSlices]]$sliceNumber - oneDat$sliceData[[1]]$sliceNumber) * fftHop
               sr <- fftLen * oneDat$maxFreq /
                   max(unlist(lapply(oneDat$sliceData, function(x) x$peakData)))
               freqRes <- sr / fftLen
           }
    )
    result$freq <- (1:ncol(bgData) + info$background[[1]]$firstBin - 1) * freqRes
    class(result) <- c('PamNoise', 'list')
    result
}

#' @export
#' @rdname loadBackgroundNoise
#' 
plotBackgroundNoise <- function(x) {
    if(!inherits(x, 'PamNoise')) {
        if(is.character(x) ||
           inherits(x, 'PamBinary')) {
            return(plotBackgroundNoise(loadBackgroundNoise(x)))
        }
        return(FALSE)
    }
    tPretty <- pretty(as.numeric(x$times), n=5)
    tLab <- convertPgDate(tPretty)
    if('freq' %in% names(x)) {
        image(x=x$times, y=x$freq, z=x$background, xlab='Time', ylab='Frequency (Hz)', xaxt='n')
        axis(1, at=tPretty, labels=tLab, cex.axis=.8)
        return(TRUE)
    }
    plot(x=x$times, y=x$background[, 1], type='l', xlab='Time', ylab='Background Level', xaxt='n')
    axis(1, at=tPretty, labels=tLab, cex.axis=.9)
    if(ncol(x$background) > 1) {
        lines(x=x$times, y=x$background[, 2], col='blue')
    }
    return(TRUE)
}