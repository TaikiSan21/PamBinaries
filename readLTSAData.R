# Reads binary data stored by the LTSA module.
#
# Inputs:
#     fid = file identifier
#     fileInfo = structure holding the file header, module header, a handle
#     data = a structure containing the standard data
# Output:
#     data = structure containing data from a single object

readLTSAData <- function(fid, fileInfo, data) {
      error <- FALSE
      a <- 127*2/log(32767)
      b <- -127
      tryCatch({
            dataLength <- pamBinRead(fid, 'int32', n=1)
            if(dataLength==0) {
                  return(list(data=data, error=error))
            }
            
            version <- fileInfo$moduleHeader$version
            
            if(version <= 1) {
                  data$startSample <- pamBinRead(fid, 'int64', n=1)
            }
            
            if(version==0) {
                  data$duration <- pamBinRead(fid, 'int64', n=1)
            }
            
            if(version <= 1) {
                  data$channelMap <- pamBinRead(fid, 'int32', n=1)
            }
            
            data$endMillis <- pamBinRead(fid, 'int64', n=1)
            data$endDate <- millisToDateNum(data$endMillis)
            data$nFFT <- pamBinRead(fid, 'int32', n=1)
            data$maxVal <- pamBinRead(fid, 'float', n=1)
            
            # Version 0 scaled the data linearly to 16 bit
            if(version==0) {
                  data$byteData <- pamBinRead(fid, 'int16', n = fileInfo$moduleHeader$fftLength / 2)
                  data$data <- data$byteData / 32767 * data$maxVal
            }
            # After version 0, the data was first scaled to 16 bit and then
            # converted to a log so that it could be saved as an 8 bit
            else {
                  data$byteData <- pamBinRead(fid, 'int8', n = fileInfo$moduleHeader$fftLength / 2)
                  data$data <- exp((data$byteData - b)/a)* data$maxVal / 32767
            }
            
            return(list(data=data, error=error))
      }, warning = function(w) {
            print(paste('Warning occurred: ', w))
            return(list(data=data, error=error))
      }, error = function(e) {
            print(paste('Error reading ', fileInfo$fileHeader$moduleType, ' data object. Data read:'))
            print(data)
            print(e)
            error <- TRUE
            return(list(data=data, error=error))
      })
}