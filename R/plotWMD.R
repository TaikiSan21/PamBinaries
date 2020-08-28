#' @title Plot Whistle Contour
#' 
#' @description Plots the entire whistle contour saved in a Pamguard Whistle &
#'   Moan Detector binary file, highlighting the selected contour
#'   
#' @param data either a \code{PamBinary} class object, or just the \code{$data} from
#'   a PamBinary object, or a single detection from the \code{$data}
#' @param id the id of the whistle to plot, either an index or Pamguard UID
#' @param \dots parameters to pass to other functions
#' 
#' @return A ggplot object
#' 
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#' 
#' @examples 
#' 
#' # load example whistle file
#' wmFile <- system.file('extdata', 'WM.pgdf', package='PamBinaries')
#' wmData <- loadPamguardBinaryFile(wmFile)
#' plotWMD(wmData, 1)
#' plotWMD(wmData, 2)
#' 
#' @importFrom dplyr distinct filter
#' @importFrom ggplot2 ggplot stat_bin_2d aes_string xlab ylab
#' @export
#' 
plotWMD <- function(data, id=1, ...) {
    if(inherits(data, 'PamBinary')) {
        return(plotWMD(data$data, id, ...))
    }
    if(is.list(data[[1]]) &&
       names(data[[1]])[1] == 'flagBitMap') {
        if(log10(id) >= 5) {
            id <- as.character(id)
        }
        return(plotWMD(data[[id]], id, ...))
    }
    if(!all(c('time', 'freq', 'contour', 'allFreq') %in% names(data))) {
        data <- contourToFreq(list(data), ...)[[1]]
    }
    # browser()
    fDiff <- unique(round(data$freq / data$contour, 3))
    if(length(fDiff) != 1) {
        stop('Inconsistent frequency resolution found, fDiff: ', fDiff)
    }
    allData <- do.call(rbind, lapply(seq_along(data$time), function(x) {
        theseFreq <- data$allFreq[, x]
        theseFreq <- theseFreq[theseFreq > 0]
        data.frame(time = data$time[x], freq = seq(from=min(theseFreq), to = max(theseFreq), by = fDiff))
    }))
    allData <- distinct(allData)
    peakData <- data.frame(time = data$time, freq = data$freq)
    tDiff <- min(unique(diff(data$time))) * .9999
    fDiff <- fDiff * .9999

    plot <- ggplot() +
        stat_bin_2d(data=allData, aes_string(x='time', y='freq'), fill='darkgrey',
                    binwidth=c(tDiff, fDiff),
                    color='black') +
        stat_bin_2d(data=peakData, aes_string(x='time', y='freq'), fill='pink',
                    binwidth=c(tDiff, fDiff),
                    color='black') +
        xlab('Time (seconds)') +
        ylab('Frequency (hertz)')
    plot
}
