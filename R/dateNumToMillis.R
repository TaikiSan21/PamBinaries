#' @title Convert Date Number to Milliseconds
#' 
#' @description Converts numeric date to millisecond date.
#'   
#' @param datenum Numeric value of a date.
#' 
#' @return Date as milliseconds
#' 
#' @author Taiki Sakai \email{taiki.sakai@noaa.gov}
#' 
#' @note Conversion to milliseconds to match how Java stores
#'   dates. Doesn't appear to ever be used.
#'   
dateNumToMillis <- function(datenum) {
    datenum * 1000
}