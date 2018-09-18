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
#' @author Taiki Sakai \email{taiki.sakai@noaa.gov}
#' 
#' @export
#'
convertPgDate <- function(dateNum) {
    as.POSIXct(dateNum, origin='1970-01-01', tz='UTC')
}