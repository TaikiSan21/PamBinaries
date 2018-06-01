#' @title Load Pamguard Binary File
#' 
#' @description This function will load in the data from a Pamguard binary file. It will
#'   figure out the type of data being read based on the header of the file.
#'   All functions based on Matlab code written by Michael Oswald.
#' 
#' @param fileName The name of the binary file to be read
#' @param \dots Arguments passed to other functions
#' @return This function returns a list containing two objects. Data contains
#'   all the binary data read. fileInfo contains metadata information for the file.
#' 
#' @author Taiki Sakai \email{taiki.sakai@noaa.gov}
#' 
#' @export
#' 
loadPamguardBinaryFile <- function(fileName, ...) {
    tryCatch({
        fid <- file(fileName, open='rb')
        
        # initialize variables
        prevPos <- -1
        dataSet <- list()
        fileInfo <- list()
        fileInfo$readModuleHeader <- readStdModuleHeader
        fileInfo$readModuleFooter <- readStdModuleFooter
        
        # main loop
        while(TRUE) {
            
            # if for some reason we're stuck at one byte, warn the user and
            # abort
            pos <- seek(fid)
            if(pos == prevPos) {
                print(paste('File stuck at byte ', pos))
                break
            }
            prevPos <- pos
            
            # read in the length of the object and the type. Move the file
            # pointer back to the beginning of the length variable, and switch
            # on the type. If we couldn't read nextLen or nextType, assume
            # that means we've hit the end of the file and break out of loop
            # browser()
            nextLen <- pamBinRead(fid, 'int32', n=1)
            nextType <- as.character(pamBinRead(fid, 'int32', n=1))
            if((length(nextLen) == 0) | (length(nextType) == 0)) {
                break
            }
            seek(fid, -8, origin='current')
            switch(nextType,
                   '-1' = {
                       fileInfo$fileHeader <- readFileHeader(fid)
                       # print(1)
                       switch(fileInfo$fileHeader$moduleType,
                              'AIS Processing' = {
                                  fileInfo$objectType <- 0
                                  fileInfo$readModuleData <- readAISData
                              },
                              'Click Detector' = {
                                  fileInfo$objectType <- 1000
                                  fileInfo$readModuleData <- readClickData
                                  fileInfo$readModuleFooter <- readClickFooter
                              },
                              'Clip Generator' = {
                                  fileInfo$objectType <- c(1, 2)
                                  fileInfo$readModuleData <- readClipData
                              },
                              'DbHt' = {
                                  fileInfo$objectType <- 1
                                  fileInfo$readModuleData <- readDbHtData
                              },
                              'DIFAR Processing' = {
                                  fileInfo$objectType <- 0
                                  fileInfo$readModuleData <- readDifarData
                              },
                              'LTSA' = {
                                  fileInfo$objectType <- 1
                                  fileInfo$readModuleHeader <- readLTSAHeader
                                  fileInfo$readModuleData <- readLTSAData
                              },
                              'Noise Monitor' = {
                                  fileInfo$objectType <- 1
                                  fileInfo$readModuleHeader <- readNoiseMonHeader
                                  fileInfo$readModuleData <- readNoiseMonData
                              },
                              'Noise Band' = {
                                  fileInfo$objectType <- 1
                                  fileInfo$readModuleHeader <- readNoiseMonHeader
                                  fileInfo$readModuleData <- readNoiseMonData
                              },
                              'NoiseBand' = {
                                  fileInfo$objectType <- 1
                                  fileInfo$readModuleData <- readNoiseBandData
                              },
                              'RW Edge Detector' = {
                                  fileInfo$objectType <- 0
                                  fileInfo$readModuleData <- readRWEDetectorData
                                  print('Right Whale Edge Detector binary file detected')
                                  print('Note that the low, high and peak frequencies are actually')
                                  print('saved as FFT slices. In order to convert values to Hz, they')
                                  print('must be multiplied by (sampleRate/fftLength)')
                              },
                              'WhistlesMoans' = {
                                  fileInfo$objectType <- 2000
                                  fileInfo$readModuleHeader <- readWMDHeader
                                  fileInfo$readModuleData <- readWMDData
                              },
                              {
                                  print(paste("Don't recognize type ", 
                                              fileInfo$fileHeader$moduleType,
                                              '. Aborting load.'))
                                  break
                              })
                   },
                   # Case 2: File Footer. The file version should have been set
                   # when we read the file header. If the file header is empty,
                   # something has gone wrong so warn the user and exit.
                   '-2' = {
                       # print(2)
                       if(length(fileInfo$fileHeader)==0){
                           print('Error: found file footer before file header. Aborting load.')
                           break
                       }
                       fileInfo$fileFooter <- readFileFooterInfo(fid, fileInfo$fileHeader$fileFormat)
                   },
                   # Case 3: Module Header. The Correct function handle should 
                   # have been set when we read the file header. If the file 
                   # header is empty, something has gone wrong so warn the user
                   # and exit
                   '-3' = {
                       # print(3)
                       if(length(fileInfo$fileHeader)==0) {
                           print('Error: found module header before file header. Aborting load.')
                           break
                       }
                       fileInfo$moduleHeader <- fileInfo$readModuleHeader(fid)
                   },
                   # Case 4: Module Footer. The correct function handle should
                   # have been set when we read the file header. If the file
                   # header is empty, something has gone wrong so warn the user
                   # and exit
                   '-4' = {
                       # print(4)
                       if(length(fileInfo$fileHeader)==0) {
                           print('Error: found module footer before file header. Aborting load.')
                           break
                       }
                       fileInfo$moduleFooter <- fileInfo$readModuleFooter(fid)
                   },
                   # Case 5: Data. The correct function handle should have been
                   # set when we read in the file header. If the file header is
                   # empty, something has gone wrong so warn the user and exit
                   {
                       # print(5)
                       if(length(fileInfo$fileHeader)==0) {
                           print('Error: found data before file header. Aborting load.')
                           break
                       }
                       dataPoint <- readPamData(fid, fileInfo, ...)
                       
                       dataSet[[length(dataSet)+1]] <- dataPoint
                   }
            )
        }
        close(fid)
        list(data=dataSet, fileInfo=fileInfo)
    }, error = function(e) {
        cat('Error reading file ', fileName)
        print(e)
    })
}

