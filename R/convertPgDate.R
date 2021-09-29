#' @title Convert Pamguard Numeric Date to POSIXct
#' 
#' @description a simple helper to convert Pamguard's numeric date to
#'   POSIXct format
#' 
#' @param dateNum date as a numeric, seconds since 1970-01-01 per standard
#'   Pamguard output. Timezone is UTC
#'   
#' @return A POSIXct date in UTC
#' 
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#' 
#' @examples 
#' 
#' # load the example click binary data, leaving date as numeric
#' clickFile <- system.file('extdata', 'Click.pgdf', package='PamBinaries')
#' clickData <- loadPamguardBinaryFile(clickFile, convertDate = FALSE)
#' # convert date to POSIXct
#' convertPgDate(clickData$data[[1]]$date)
#' 
#' @export
#'
convertPgDate <- function(dateNum) {
    if(is.null(dateNum) ||
       all(is.na(dateNum))) {
        return(dateNum)
    }
    if(!is.numeric(dateNum)) {
        stop('Date must be numeric, "', dateNum, '" is class ', class(dateNum))
    }
    if(any(is.na(dateNum))) {
        dateNum[!is.na(dateNum)] <- as.POSIXct(dateNum[!is.na(dateNum)], origin='1970-01-01', tz='UTC')
        return(dateNum)
    }
    as.POSIXct(dateNum, origin='1970-01-01', tz='UTC')
}
