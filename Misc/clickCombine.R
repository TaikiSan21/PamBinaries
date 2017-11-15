library(data.table)
library(dplyr)
library(PamBinaries)

###################################################################################################
# Can either specify specific folders to process, or a directory containing all folders to process.
# Will write files into the current directory. The default settings are to skip reading wave
# files (getWave=FALSE), process all folders in the current directory, and print progress
# updates every 10 files. It will also create a ClickErrorLog.txt in the current directory if
# any errors are encountered that shows what the errors were. It should (hopefully) keep running even
# if it encounters an error reading a file, then you can use the log to find out which files you need
# to go back and run again.
###################################################################################################

# EXAMPLES: 
# To create a CSV for each folder in directory 'S:/AllMyBinaries/' :
# writeClickFolder(dir='S:/AllMyBinaries/', getWave=FALSE)
#
# To read only specific folders c('S:/AllMyBinaries/Folder1', 'S:/AllMyBinaries/Folder2') 
# (Useful if you only need to run a subset of folders in a directory) :
# writeClickFolder(folders=c('S:/AllMyBinaries/Folder1', 'S:/AllMyBinaries/Folder2'))
#
# You can set the frequency of message updates (ie. 'On file 1 of 100') by setting
# messageFrequency, or tell it to be quiet with quiet=TRUE:
#
# Only print progress message every 20 files:
# writeClickFolder(dir='S:/AllMyBinaries/', messageFrequency=20)
#
# Or no messages:
# writeClickFolder(dir='S:/AllMyBinaries/', quiet=TRUE)

writeClickFolder <- function(dir='.', folders=NULL, ...) {
    if(is.null(folders)) {
        foldersToCombine <- list.dirs(dir, recursive=FALSE)
    } else {
        foldersToCombine <- folders
    }
    lapply(
        foldersToCombine, function(f) {
            AllClicks <- combineClickFiles(f, ...)
            folderName <- gsub('.*(/|\\\\)', '', f)
            write.csv(rbindlist(lapply(AllClicks, function(x) {
                x$clickDf})), file=paste0(folderName, '.CSV'),
                row.names=FALSE)
        })
}

combineClickFiles <- function(folder, fileList=NULL, quiet=FALSE, getWave=FALSE, messageFrequency=10, clickNos=NULL) {
    hSens <- -172.1
    gain <- 1
    adcPeakPeak <- 2.828
    if(is.null(fileList)) {
        files <- list.files(folder, pattern='Click_Detector_Click_Detector_Clicks.*\\.pgdf$',
                            recursive=TRUE, full.names=TRUE)
    } else {
        files <- fileList
    }
    
    logName <- 'ClickErrorLog.txt'
    nFiles <- length(files)
    iFile <- 1
    errors <- c()
    
    AllClicksList <- lapply(files,
                            function(file) {
                                # Show progress
                                result <- list()
                                shortFile <- sub('.*(([/|\\\\].*){3})$', '\\1', file)
                                if(!quiet) {
                                    # print(mem_used())
                                    # time <- proc.time()
                                    if((iFile %% messageFrequency) == 1) {
                                        cat(iFile, ' of', nFiles, '  ', shortFile, '\n')
                                    }
                                }
                                iFile <<- iFile + 1
                                # Gather clicks
                                tryCatch({
                                    # log <- file(logName, open = 'wt')
                                    # sink(log, type='message')
                                    binaryData <- loadPamguardBinaryFile(file, getWave=getWave)
                                    #browser()
                                    nClicks <- length(binaryData$data)
                                    if(nClicks > 0) {
                                        if(!getWave) {
                                            result$clickDf <- tryCatch({
                                                rbindlist(lapply(binaryData$data, function(d) {
                                                    class(d) <- 'data.frame'
                                                    attr(d, 'row.names') <- .set_row_names(length(d[[1]]))
                                                    d}))}, error=function(e) {
                                                        print('Fast version failed, going slow.')
                                                        do.call(rbind, lapply(binaryData$data, function(d) {
                                                            data.frame(d)
                                                        }
                                                        ))}) %>% mutate(ClickNo = 0:(n()-1),
                                                                        BinaryFile = shortFile)
                                            # Gather waves
                                        } else {
                                            waves <- vector('list', nClicks)
                                            for(i in seq_along(waves)) {
                                                waves[[i]] <- binaryData$data[[i]]$wave
                                            }
                                            result <- waves
                                            names(result) <- 0:(nClicks-1)
                                        }
                                        # if(iFile==2) 1+'1'
                                        
                                    }
                                    # if(!quiet) {
                                    #     cat('\n Took ', (proc.time()-time)[3],'\n')
                                    # }
                                    if(!is.null(clickNos)) {
                                        result <- result[(clickNos+1)]
                                    }
                                    result
                                }, error = function(e) {
                                    errors <<- c(errors, shortFile)
                                    cat(paste('Error in file ', shortFile, ': \n', e), file=logName, append=TRUE)
                                    list(clickDf=data.frame(), wave=list())
                                })
                            }
                            
    )
    fileNames <- sub('.*([/|\\\\](.*))$', '\\2', files)
    if(length(errors) > 0) {
        cat('WARNING: Encountered ', length(errors), ' error(s). Unable to load files: \n',
            errors, '\n see log file ', logName, ' for details.')
        # errorFiles <- sub('.*([/|\\\\](.*))$', '\\2', errors)
        # fileNames <- fileNames[!(fileNames %in% errorFiles)]
    }
    names(AllClicksList) <- fileNames
    AllClicksList
}
