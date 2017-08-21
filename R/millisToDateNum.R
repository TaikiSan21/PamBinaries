#' @title Convert Java Millisecond Time to R
#' 
#' @description Converts Java millisecond time into numeric
#'   time that R uses. 
#'   
#' @param millis Millisecond time from Java
#' 
#' @return Numeric time used by R.
#' 
#' @author Taiki Sakai \email{taiki.sakai@noaa.gov}
#' 
#' @note Original function was more relevant as Matlab and Java
#'   use different time origins. Java & R both use 1970-01-01,
#'   but Java stores as milliseconds vs seconds in R.
#'   
millisToDateNum <- function(millis) {
    millis / 1000
}