library(data.table)
library(dplyr)
# library(PamBinaries)
# setwd('./TestFiles')

taikiClicks <- function(dir='.', folders=NULL, outDir=NULL, ...) {
    if(is.null(folders)) {
        foldersToCombine <- list.dirs(dir, recursive=FALSE)
    } else {
        foldersToCombine <- folders
    }
    if(is.null(outDir)) {
        outDir <- dir
    }
    lapply(
        foldersToCombine, function(f) {
            AllClicks <- combineClickFiles(f, ...)
            folderName <- gsub('.*(/|\\\\)', '', f)
            write.csv(rbindlist(lapply(AllClicks, function(x) {
                x$clickDf})), file=paste0('TaikiComp',folderName, '.CSV'),
                row.names=FALSE)
        })
}

combineClickFiles <- function(folder, fileList=NULL, quiet=FALSE, getWave=FALSE, messageFrequency=10) {
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
                                    log <- file(logName, open = 'wt')
                                    sink(log, type='message')
                                binaryData <- loadPamguardBinaryFile(file, getWave=getWave)
                                #browser()
                                nClicks <- length(binaryData$data)
                                if(nClicks > 0) {
                                    if(!getWave) {
                                        result$clickDf <- rbindlist(lapply(binaryData$data, function(d) {
                                            tryCatch({
                                                class(d) <- 'data.frame'
                                                attr(d, 'row.names') <- .set_row_names(length(d[[1]]))
                                                d}, error=function(e) {
                                                    print('Fast version failed, going slow.')
                                                    data.frame(d)
                                                })
                                        })) %>% mutate(ClickNo = 0:(n()-1),
                                                       BinaryFile = shortFile)
                                        # Gather waves
                                    } else {
                                        waves <- vector('list', nClicks)
                                        for(i in seq_along(waves)) {
                                            waves[[i]] <- binaryData$data[[i]]$wave
                                        }
                                        result$wave <- waves
                                    }
                                    # if(iFile==2) 1+'1'
                                    
                                }
                                # if(!quiet) {
                                #     cat('\n Took ', (proc.time()-time)[3],'\n')
                                # }
                                result
                                }, error = function(e) {
                                    errors <<- c(errors, shortFile)
                                    cat(paste('Error in file ', shortFile, ': \n', e), file=logName, append=TRUE)
                                    list(clickDf=data.frame(), wave=matrix())
                                })
                            }
                            
    )
    if(length(errors) > 0) {
        cat('WARNING: Encountered ', length(errors), ' error(s). Unable to load files: \n',
            errors, '\n see log file ', logName, ' for details.')
    }
    
    AllClicksList
}


