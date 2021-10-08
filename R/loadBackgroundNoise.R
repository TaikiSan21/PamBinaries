#' @title Load and Format Background Noise Data
#' 
#' @description Reads and formats background noise data from Pamguard binary files
#'   or if not present in the original file will try to read the accompanying 
#'   .pgnf noise file if it exists
#' 
#' @param x character pointing to a Pamguard binary file, or a \code{PamBinary}
#'   object created by \link{loadPamguardBinaryFile}. For plotting or combining, 
#'   either of these or the output from \link{loadBackgroundNoise}
#' @param forPlot logical flag when combining noise data. If used for plotting
#'   purposes this will insert NA columns into background data so that images show
#'   up with time gaps as expected. Leave as \code{FALSE} unless you are sure you
#'   want this.
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
#' @importFrom stats sd
#' @export
#'
loadBackgroundNoise <- function(x) {
    if(is.character(x) &&
       file.exists(x)) {
        x <- loadPamguardBinaryFile(x, skipLarge=TRUE)
    }
    info <- x$fileInfo
    oneDat <- x$data[[1]]
    if(is.null(info$background)) {
        noiseFile <- gsub('df$', 'nf', info$fileName)
        if(file.exists(noiseFile)) {
            noise <- loadPamguardBinaryFile(noiseFile)
            noise$data <- x$data
            if(is.null(noise$fileInfo$background)) {
                return(NULL)
            }
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
    result <- list(detector=gsub(' ', '_', info$fileHeader$moduleName),
                   times = times, background=bgData)
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
#' @importFrom graphics title
#' @rdname loadBackgroundNoise
#' 
plotBackgroundNoise <- function(x) {
    x <- combineBackgroundNoise(x, forPlot=TRUE)
    for(i in seq_along(x)) {
        tPretty <- pretty(as.numeric(x[[i]]$times), n=5)
        tLab <- convertPgDate(tPretty)
        if('freq' %in% names(x[[i]])) {
            # Adjusting background data to plot more clearly
            bgFix <- x[[i]]$background
            noPlot <- 1:2
            # if(grepl('[Cc]epstrum', names(x)[i])) {
            #     bgFix[, 1] <- max(bgFix[, -1], na.rm=TRUE)
            # }
            bgFix[, noPlot] <- max(bgFix[, -noPlot], na.rm=TRUE)
            lim <- mean(bgFix, na.rm=TRUE) + c(-1,1) * 3 * sd(bgFix, na.rm=TRUE)
            bgFix[bgFix < lim[1]] <- lim[1]
            bgFix[bgFix > lim[2]] <- lim[2]
            image(x=x[[i]]$times, y=x[[i]]$freq, z=bgFix, xlab='Time', ylab='Frequency (Hz)', xaxt='n')
            axis(1, at=tPretty, labels=tLab, cex.axis=.8)
            title(main=paste0(names(x)[i], ' Background Noise'))
            next
        }
        plot(x=x[[i]]$times, y=x[[i]]$background[, 1], type='l', xlab='Time', ylab='Background Level', xaxt='n')
        axis(1, at=tPretty, labels=tLab, cex.axis=.9)
        title(main=paste0(names(x)[i], ' Background Noise'))
        if(ncol(x[[i]]$background) > 1) {
            lines(x=x[[i]]$times, y=x[[i]]$background[, 2], col='blue')
        }
    }
    invisible(TRUE)
}

#' @export
#' @rdname loadBackgroundNoise
#' 
combineBackgroundNoise <- function(x, forPlot=FALSE) {
    if(inherits(x, 'PamNoise')) {
        x <- list(x)
    }
    if(all(is.character(x))) {
        x <- lapply(x, loadBackgroundNoise)
    }
    x <- x[!sapply(x, is.null)]
    for(i in seq_along(x)) {
        names(x)[i] <- x[[i]]$detector
        keep <- x[[i]][2:length(x[[i]])]
        x[[i]] <- keep
    }
    detNames <- unique(names(x))
    result <- vector('list', length = length(detNames))
    names(result) <- detNames
    for(n in detNames) {
        whichThisName <- which(names(x) == n)
        thisData <- x[whichThisName]
        if(forPlot &&
           length(thisData) > 1) {
            for(i in seq_along(thisData)) {
                if(i == 1) {
                    thisData[[i]]$background <- rbind(thisData[[i]]$background, NA)
                    thisData[[i]]$times <- c(thisData[[i]]$times, thisData[[i]]$times[length(thisData[[i]]$times)]+.01)
                } else {
                    thisData[[i]]$background <- rbind(NA, thisData[[i]]$background)
                    thisData[[i]]$times <- c(thisData[[i]]$times[1]-.01, thisData[[i]]$times)
                }
            }
        }
        result[[n]] <- list(times = unlist(lapply(thisData, function(w) w$times), use.names=FALSE),
                            background = do.call(rbind, lapply(thisData, function(w) w$background)))
        if('freq' %in% names(thisData[[1]])) {
            result[[n]]$freq <- thisData[[1]]$freq
        }
    }
    # x <- squishList(x)
    for(i in seq_along(result)) {
        if('freq' %in% names(result[[i]])) {
            result[[i]]$freq <- result[[i]]$freq[1:ncol(result[[i]]$background)]
        }
        if(!inherits(result[[i]]$times, 'POSIXct')) {
            result[[i]]$times <- convertPgDate(result[[i]]$times)
        }
        dupeTime <- duplicated(result[[i]]$times)
        result[[i]]$times <- result[[i]]$times[!dupeTime]
        timeSort <- sort(as.numeric(result[[i]]$times), index.return=TRUE)$ix
        result[[i]]$times <- result[[i]]$times[timeSort]
        result[[i]]$background <- result[[i]]$background[!dupeTime, ][timeSort, ]
        class(result[[i]]) <- c('PamNoise', 'list')
    }
    result
}
