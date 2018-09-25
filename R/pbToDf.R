#' @title Convert a PamBinary Object to Data Frame
#' 
#' @description Converts a PamBinary object into a data frame. The data.frame
#'   will combine all of the data from the \code{data} part of the PamBinary
#'   object, but will not include annotations data, click waveforms, DIFAR demux
#'   data, or contours from the WMD detector. These are skipped because they are
#'   either inconsistent in their size, or are large objects. The function 
#'   \code{pbToDf} is also called when \code{as.data.frame} is called on a PamBinary
#'   class object.
#'   
#' @param pb a PamBinary class object created by 
#'   \code{\link[PamBinaries]{loadPamguardBinaryFile}}
#' @param \dots Unused parameters passed in from other methods
#' 
#' @return a data.frame containing most of the binary data read in. Will not
#'   contain annotation data, click waveforms, DIFAR demux data, or contour
#'   information from WMD detector. These are skipped because they are either
#'   incosistent in their size, or are large objects.
#' 
#' @author Taiki Sakai \email{taiki.sakai@noaa.gov}
#' 
#' @export
#' 
pbToDf <- function(pb) {
    skip <- c('annotations', 'wave', 'contour', 'contWidth', 'sliceData', 'demuxData')
    # Case 1: either a PamBinary class, or has the right pieces but not the class
    if('PamBinary' %in% class(pb) ||
       all(c('data', 'fileInfo') %in% names(pb))) {
        return(
            do.call(rbind, lapply(pb$data, function(x) {
                data.frame(x[!(names(x) %in% skip)])
            }))
        )
    }
    # Case 2: Just the $data part, we should do it anyway but warn.
    if(all(c('flagBitMap', 'identifier') %in% names(pb[[1]]))) {
        warning('It looks like you input just the $data, please use',
                'the entire "PamBinary" output next time.')
        return(
            do.call(rbind, lapply(pb, function(x) {
                data.frame(x[!(names(x) %in% skip)])
            }))
        )
    }
    stop("Input doesn't look like a PamBinary output.")
}

#' @export
#' 
as.data.frame.PamBinary <- function(x, ...) {
    pbToDf(x)
}
