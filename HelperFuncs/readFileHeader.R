# Read file headers

readFileHeader <- function(file, readExtra=FALSE) {
    header <- list()
    # browser()
    header$length <- pamBinRead(file, 'int32', n=1)
    header$identifier <- pamBinRead(file, 'int32', n=1)
    header$fileFormat <- pamBinRead(file, 'int32', n=1)
    header$pamguard <- pamBinRead(file, 'character', n=12)
    header$version <- readJavaUTFString(file)$str # These arent marked as only giving string part in matlab??
    header$branch <- readJavaUTFString(file)$str
    ######################  Need to look into R's specification of dates
    header$dataDate <- millisToDateNum(pamBinRead(file, 'int64', n=1))
    header$analysisDate <- millisToDateNum(pamBinRead(file, 'int64', n=1))
    ###################
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
