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
#' @param templateNames if using the click template classifier, the names of the species
#'   for the click templates. These will be used as the names of the columns in the
#'   dataframe, and the length of this must exactly match the number of templates used
#' @param \dots Unused parameters passed in from other methods
#' 
#' @return a data.frame containing most of the binary data read in. Will not
#'   contain annotation data, click waveforms, DIFAR demux data, or contour
#'   information from WMD detector. These are skipped because they are either
#'   incosistent in their size, or are large objects.
#' 
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#' 
#' @importFrom dplyr bind_rows
#' @export
#' 
pbToDf <- function(pb, templateNames = NULL) {
    skip <- c('annotations', 'wave', 'contour', 'contWidth', 'sliceData', 'demuxData')
    good <- FALSE
    # Case 1: either a PamBinary class, or has the right pieces but not the class
    if('PamBinary' %in% class(pb) ||
       all(c('data', 'fileInfo') %in% names(pb))) {
        justData <- pb$data
        good <- TRUE
    }
    # Case 2: Just the $data part, we should do it anyway but warn.
    if(all(c('flagBitMap', 'identifier') %in% names(pb[[1]]))) {
        warning('It looks like you input just the $data, please use',
                'the entire "PamBinary" output next time.')
        justData <- pb
        good <- TRUE
    }
    keepIx <- !(names(justData[[1]]) %in% skip)
    if(good) {
        return(
            bind_rows(lapply(justData, function(x) {
                if(is.null(templateNames)) {
                    x[keepIx]
                } else {
                    ct <- x$annotations$mclassification$threshold
                    names(ct) <- templateNames
                    c(x[keepIx], ct)
                }
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
