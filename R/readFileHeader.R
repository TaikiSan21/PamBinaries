#' @title Read File Header
#' 
#' @description Reads file header information common to all files
#'   
#' @param file binary file to be read
#' @param readExtra flag if there is extra information to read
#' 
#' @return header information common to all files
#' 
#' @author Taiki Sakai \email{taiki.sakai@noaa.gov}
#' 
readFileHeader <- function(file, readExtra=FALSE) {
    header <- list()
    header$length <- pamBinRead(file, 'int32', n=1)
    header$identifier <- pamBinRead(file, 'int32', n=1)
    header$fileFormat <- pamBinRead(file, 'int32', n=1)
    header$pamguard <- pamBinRead(file, 'character', n=12)
    header$version <- readJavaUTFString(file)$str
    header$branch <- readJavaUTFString(file)$str
    header$dataDate <- millisToDateNum(pamBinRead(file, 'int64', n=1))
    header$analysisDate <- millisToDateNum(pamBinRead(file, 'int64', n=1))
    header$startSample <- pamBinRead(file, 'int64', n=1)
    header$moduleType <- readJavaUTFString(file)$str
    header$moduleName <- readJavaUTFString(file)$str
    header$streamName <- readJavaUTFString(file)$str
    header$extraInfoLen <- pamBinRead(file, 'int32', n=1)
    if(readExtra){
        header$extraInfo <- pamBinRead(file, 'int8', n=header$extraInfoLen)
    } else {
        seek(file, header$extraInfoLen, origin='current')
    }
    return(header)
}
