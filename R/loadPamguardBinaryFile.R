#' @title Load Pamguard Binary File
#' 
#' @description This function will load in the data from a Pamguard binary file. It will
#'   figure out the type of data being read based on the header of the file.
#'   All functions based on Matlab code written by Michael Oswald.
#' 
#' @param fileName The name of the binary file to be read
#' @param skipLarge Should we skip large parts of binaries? Currently only applicable
#'   to whistle, click, and DIFAR data
#' @param skipData Should we skip all data and only read headers and footers?
#' @param debug logical flag to show more info on errors
#' @param keepUIDs If not \code{NULL}, a vector of UIDs to read. All UIDs not in this
#'   vector will not be read.
#' @param convertDate logical flag to convert date from numeric to POSIXct. Defaults to 
#'   \code{FALSE} for speed, can reduce time by 
#' @param \dots Arguments passed to other functions
#' @return This function returns a list containing two objects. Data contains
#'   all the binary data read. fileInfo contains metadata information for the file.
#' 
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#' 
#' @examples 
#' 
#' # read example whistle data
#' wmFile <- system.file('extdata', 'WM.pgdf', package='PamBinaries')
#' whistleData <- loadPamguardBinaryFile(wmFile)
#' # works the same for different kinds of binary files
#' clickFile <- system.file('extdata', 'Click.pgdf', package='PamBinaries')
#' clickData <- loadPamguardBinaryFile(clickFile)
#' # convert date to POSIXct (default does not because it is faster)
#' clickPOSIX <- loadPamguardBinaryFile(clickFile, convertDate = TRUE)
#' clickData$data[[1]]$date
#' clickPOSIX$data[[1]]$date
#' # read only the fileInfo portion, has empty $data item
#' clickInfo <- loadPamguardBinaryFile(clickFile, skipData = TRUE)
#' # skip reading the large click waveforms, much faster if you dont need them
#' clickLess <- loadPamguardBinaryFile(clickFile, skipLarge = TRUE)
#' object.size(clickData)
#' object.size(clickLess)
#' # only read specific UID numbers
#' clickSpecific <- loadPamguardBinaryFile(clickFile, keepUIDs = c(4000006, 4000007))
#' names(clickSpecific$data)
#' 
#' @export
#' 
loadPamguardBinaryFile <- function(fileName, skipLarge=FALSE, skipData=FALSE, 
                                   debug=FALSE, keepUIDs=NULL, convertDate=FALSE, ...) {
    tryCatch({
        fid <- file(fileName, open='rb')
        
        # initialize variables
        nBackground <- 0
        prevPos <- -1
        dataSet <- list()
        fileInfo <- list()
        backgroundData <- list()
        fileInfo$fileName <- fileName
        fileInfo$readModuleHeader <- readStdModuleHeader
        fileInfo$readModuleFooter <- readStdModuleFooter
        doneUIDs <- character(0)
        
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
                                  switch(fileInfo$fileHeader$streamName,
                                         'Clicks' = {
                                             fileInfo$objectType <- 1000
                                             fileInfo$readModuleData <- readClickData
                                             fileInfo$readModuleFooter <- readClickFooter
                                             fileInfo$readBackgroundData <- readClickBackground
                                         },
                                         'Trigger Background' = {
                                             fileInfo$objectType <- 0
                                             fileInfo$readModuleData <- readClickTriggerData
                                             fileInfo$readModuleHeader <- readClickTriggerHeader
                                         })
                              },
                              'SoundTrap Click Detector' = {
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
                              'Deep Learning Classifier' = {
                                  switch(fileInfo$fileHeader$streamName,
                                         'DL_detection' = {
                                             fileInfo$objectType <- 1
                                             fileInfo$readModuleData <- readDLDetData
                                         },
                                         'DL_Model_Data' = {
                                             # TODO
                                         },
                                         'DL Model Data' = {
                                             # TODO same above
                                         }
                                  )
                              },
                              'GPL Detector' = {
                                  switch(fileInfo$fileHeader$streamName,
                                         'GPL Detections' = {
                                             fileInfo$readModuleData <- readGPLDetections
                                             fileInfo$readBackgroundData <- readSpectralBackground
                                         })
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
                              'Gemini Threshold Detector' = {
                                  fileInfo$objectType <- 0
                                  fileInfo$readModuleData <- readTritechTrack
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
                                  fileInfo$readBackgroundData <- readSpectralBackground
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
                       if(skipData) {
                           seek(fid, pos + nextLen, origin='start')
                           next
                       }
                       if(length(fileInfo$fileHeader)==0) {
                           print('Error: found data before file header. Aborting load.')
                           break
                       }
                       dataPoint <- readPamData(fid=fid, fileInfo=fileInfo, skipLarge=skipLarge, 
                                                debug=debug, keepUIDs=keepUIDs, ...)
                       if(!is.null(dataPoint)) {
                           if (nextType == -6) {
                               nBackground <- nBackground + 1
                               backgroundData[[nBackground]] <- dataPoint
                               next
                           }
                           if(dataPoint$UID %in% doneUIDs) {
                               next
                           }
                           dataSet[[length(dataSet)+1]] <- dataPoint
                           doneUIDs <- c(doneUIDs, dataPoint$UID)
                       }
                       # Stop if at end of list of UIDs
                       # MAYBE PROBLEM: This skips the footers.
                       if(length(keepUIDs) > 0 &&
                          all(keepUIDs %in% doneUIDs)) {
                           skipData <- TRUE
                       }
                   }
            )
        }
        close(fid)
        
        if(length(dataSet) > 0 && 
            'UID' %in% names(dataSet[[1]])) {
            names(dataSet) <- sapply(dataSet, function(x) x$UID)
        }
        if(convertDate) {
            for(i in seq_along(dataSet)) {
                dataSet[[i]]$date <- convertPgDate(dataSet[[i]]$date)
            }
        }
        if (nBackground > 0) {
            fileInfo$background <- backgroundData
        }
        
        result <- list(data=dataSet, fileInfo=fileInfo)
        class(result) <- c('PamBinary', 'list')
        result
    }, error = function(e) {
        cat('Error reading file ', fileName)
        print(e)
        NULL
    })
}

