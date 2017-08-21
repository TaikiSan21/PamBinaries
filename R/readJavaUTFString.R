#' @title Read Java UTF-8 String
#' 
#' @description Reads a Java UTF-8 string. The first 2 bytes are the
#'   length of the string, then the string itself.
#'   
#' @param file binary file to be read
#' 
#' @return the string and its length
#' 
#' @author Taiki Sakai \email{taiki.sakai@noaa.gov}
#' 
readJavaUTFString <- function(file) {
    len <- pamBinRead(file, 'int16', n=1)
    str <- pamBinRead(file, 'character', n=len)
    list(len=len, str=str)
}