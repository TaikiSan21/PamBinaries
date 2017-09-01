#' @title Read Pamguard Data
#' 
#' @description Reads in the object data that is common to all modules. This 
#'   reads up to (but not including) the object binary length, and then calls 
#'   a function to read the module-specific data.
#'   
#' @param fid binary file identifier
#' @param fileInfo structure holding the file header, module header, and the
#'   appropriate function to read module specific data
#' 
#' @return a structure containing data from a single object
#' 
#' @author Taiki Sakai \email{taiki.sakai@noaa.gov}
#' 
readPamData <- function(fid, fileInfo, ...) {
    ### UNSURE OF WHAT THE RESULTS ARE IN CASE OF ERROR ###
    # set constants to match flag bitmap constants in class
    # DataUnitBaseData.java. The following constants match header version 4.
    TIMEMILLIS        = strtoi('1', base=16)
    TIMENANOS         = strtoi('2', base=16)
    CHANNELMAP        = strtoi('4', base=16)
    UID               = strtoi('8', base=16)
    STARTSAMPLE       = strtoi('10', base=16)
    SAMPLEDURATION    = strtoi('20', base=16)
    FREQUENCYLIMITS   = strtoi('40', base=16)
    MILLISDURATION    = strtoi('80', base=16)
    TIMEDELAYSSECS    = strtoi('100', base=16)
    
    # initialize a new variable to hold the data
    data <- list()
    data$flagBitMap <- 0
    
    # caclulate where the next object starts, in case there is an error trying
    # to read this one
    objectLen <- pamBinRead(fid, 'int32', n=1)
    curObj <- seek(fid)
    nextObj <- curObj + objectLen
    
    # first thing to check is that this is really the type of object we think
    # it should be, based on the file header. If not, warn the user, move the
    # pointer to the next object, and exit
    data$identifier <- pamBinRead(fid, 'int32', n=1)
    # browser()
    if(any(data$identifier == fileInfo$objectType)) {
        # Do nothing here- couldn't figure out a clean way of checking if
        # number wasn't in array
    } else {
        print(paste('Error - Object Identifier does not match ',
                    fileInfo$fileHeader$moduleType,
                    ' type. Aborting data read.'))
        seek(fid, nextObj, origin='start')
        return(list(fid=fid, fileInfo=fileInfo))
    }
    
    # Read the data, starting with the standard data that every data unit has
    version <- fileInfo$fileHeader$fileFormat
    tryCatch({
        data$millis <- pamBinRead(fid, 'int64', n=1)
        
        if(version >= 3) {
            data$flagBitmap <- pamBinRead(fid, 'int16', n=1)
        }
        
        if((version == 2) | (bitwAnd(data$flagBitMap, TIMENANOS) != 0)) {
            data$timeNanos <- pamBinRead(fid, 'int64', n=1)
        }
        
        if((version==2) | (bitwAnd(data$flagBitMap, CHANNELMAP) != 0)) {
            data$channelMap <- pamBinRead(fid, 'int32', n=1)
        }
        
        if(bitwAnd(data$flagBitMap, UID)==UID) {
            data$UID <- pamBinRead(fid, 'int64', n=1)
        }
        
        if(bitwAnd(data$flagBitMap, STARTSAMPLE) != 0) {
            data$startSample <- pamBinRead(fid, 'int64', n=1)
        }
        
        if(bitwAnd(data$flagBitMap, SAMPLEDURATION) != 0) {
            data$sampleDuration <- pamBinRead(fid, 'int32', n=1)
        }
        
        if(bitwAnd(data$flagBitMap, FREQUENCYLIMITS) != 0) {
            minFreq <- pamBinRead(fid, 'float', n=1)
            maxFreq <- pamBinRead(fid, 'float', n=1)
            data$freqLimits <- c(minFreq, maxFreq)
        }
        
        if(bitwAnd(data$flagBitMap, MILLISDURATION) != 0) {
            data$millisDuration <- pamBinRead(fid, 'float', n=1)
        }
        
        if(bitwAnd(data$flagBitMap, TIMEDELAYSSECS) != 0) {
            data$numTimeDelays <- pamBinRead(fid, 'int16', n=1)
            td <- rep(0, data$numTimeDelays)
            for(i in 1:data$numTimeDelays) {
                td[i] <- pamBinRead(fid, 'float', n=1)
            }
            data$timeDelays <- td
        }
        
        # set date, to maintain backwards compatibility
        data$date <- as.POSIXct(millisToDateNum(data$millis), origin='1970-01-01', tz='UTC')
        # data$date <- millisToDateNum(data$millis)
        # now read the module-specific data
        if(class(fileInfo$readModuleData)=='function') {
            result <- fileInfo$readModuleData(fid, fileInfo, data, ...)
            data <- result$data
            if(result$error) {
                print(paste('Error - cannot retrieve ', 
                            fileInfo$fileHeader$moduleType,
                            ' data properly.'))
                seek(fid, nextObj, origin='start')
                return(data)
            }
        }
        return(data)
    # }, warning = function(w) {
    #     print(paste('Warning occurred: ', w))
    #     return(data)
    }, error = function(e) {
        print('Error loading object data')
        print(data)
        print(e)
        seek(fid, nextObj, origin='start')
    })
}

