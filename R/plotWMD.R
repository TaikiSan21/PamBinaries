#' @title Plot Whistle Contour
#' 
#' @description Plots the entire whistle contour saved in a Pamguard Whistle &
#'   Moan Detector binary file, highlighting the selected contour
#'   
#' @param data either a \code{PamBinary} class object, or just the \code{$data} from
#'   a PamBinary object, or a single detection from the \code{$data}
#' 
#' @return A ggplot object
#' 
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#' 
#' @importFrom dplyr distinct filter
#' @importFrom ggplot2 ggplot stat_bin_2d aes
#' @export
#' 
plotWMD <- function(data, id=1) {
    if(inherits(data, 'PamBinary')) {
        return(plotWMD(data$data, id))
    }
    if(is.list(data[[1]]) &&
       names(data[[1]])[1] == 'flagBitMap') {
        if(log10(id) >= 5) {
            id <- as.character(id)
        }
        return(plotWMD(data[[id]], id))
    }
    if(!all(c('time', 'freq', 'contour', 'allFreq') %in% names(data))) {
        data <- contourToFreq(list(data))[[1]]
    }
    # browser()
    fDiff <- unique(data$freq / data$contour)
    allData <- do.call(rbind, lapply(seq_along(data$time), function(x) {
        theseFreq <- data$allFreq[, x]
        theseFreq <- theseFreq[theseFreq > 0]
        data.frame(time = data$time[x], freq = seq(from=min(theseFreq), to = max(theseFreq), by = fDiff))
    }))
    allData <- distinct(allData)
    peakData <- data.frame(time = data$time, freq = data$freq)
    tDiff <- min(unique(diff(data$time))) * .9999
    fDiff <- unique(data$freq / data$contour) * .9999

    plot <- ggplot() +
        stat_bin_2d(data=allData, aes(x=time, y=freq), fill='darkgrey',
                    binwidth=c(tDiff, fDiff),
                    color='black') +
        stat_bin_2d(data=peakData, aes(x=time, y=freq), fill='pink',
                    binwidth=c(tDiff, fDiff),
                    color='black')
    plot
}
    
        
