# Reads binary data stored by the dBHt module
#
# Inputs:
#     fid = file identifier
#     fileInfo = structure holding the file header, module header, a handle
#     data = a structure containing the standard data
# Output:
#     data = structure containing data from a single object

readDbHtData <- function(fid, fileInfo, data) {
      error <- FALSE
      
      tryCatch({
            dataLength <- pamBinRead(fid, 'int32', n=1)
            if(dataLength==0) {
                  return(list(data=data, error=error))
            }
            
            version <- fileInfo$moduleHeader$version
            
            if(version <= 1) {
                  data$startSample <- pamBinRead(fid, 'int64', n=1)
                  data$channelMap <- pamBinRead(fid, 'int32', n=1)
            }
            
            data$rms <- pamBinRead(fid, 'int16', n=1) / 100
            data$zeroPeak <- pamBinRead(fid, 'int16', n=1) / 100
            data$peakPeak <- pamBinRead(fid, 'int16', n=1) / 100
            
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