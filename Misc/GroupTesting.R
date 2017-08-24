# Comparison to Jay's
safeDf <- function(data) {
    do.call(rbind, lapply(data, function(x) {
        tooBig <- which(sapply(x, function(y) !is.null(dim(y))))
        data.frame(x[-tooBig])
    }))}

fastDf <- function(data) {
    do.call(rbind, lapply(data, function(x) {
        # tooBig <- which(sapply(x, function(y) !is.null(dim(y))))
        tooBig <- grep('wave', names(x))
        # tooBig <- 15
        x <- x[-tooBig]
        tryCatch({
            class(x) <- 'data.frame'
            attr(x, 'row.names') <- .set_row_names(length(x[[1]]))
            x}, error = function(e) {
                print(e)
                data.frame(x)
            })
    }))}
fastDt <- function(data) {
    rbindlist(lapply(data, function(x) {
        # tooBig <- which(sapply(x, function(y) !is.null(dim(y))))
        # tooBig <- 15
        tooBig <- grep('wave', names(x))
        x <- x[-tooBig]
        tryCatch({
            class(x) <- 'data.frame'
            attr(x, 'row.names') <- .set_row_names(length(x[[1]]))
            x}, error = function(e) {
                print(e)
                data.frame(x)
            })
    }))}
# setwd("S:\\1650_PASCAL_2016\\Data\\DASBR_Pamguard_Post_Processing\\Binaries\\Final_Dataset")
setwd('./TestFiles')
taikiTime <- function(fileFun, foldFun) {
    lapply(
        list.dirs(recursive=FALSE), function(folder) {
            files <- list.files(folder, pattern='Click_Detector_Click_Detector_Clicks.*\\.pgdf$',
                                recursive=TRUE, full.names=TRUE)
            nFiles <- length(files)
            iFile <- 1
            AllClicksList <- lapply(files,
                                    function(file) {
                                        # Show progress
                                        # cat(iFile, ' of', nFiles, '  ', file, '\n')
                                        iFile <<- iFile + 1
                                        
                                        # Gather clicks
                                        binaryData <- loadPamguardBinaryFile(file)
                                        clickDf <- fileFun(binaryData$data) %>%
                                            mutate(ClickNo = 0:(n()-1),
                                                   BinaryFile = sub('./', '/', file))
                                        
                                        # Gather waves
                                        nClicks <- length(binaryData$data)
                                        if(nClicks > 0) {
                                            waves <- vector('list', nClicks)
                                            for(i in seq_along(waves)) {
                                                waves[[i]] <- binaryData$data[[i]]$wave
                                            }
                                        }
                                        
                                        clickDf
                                    })
            AllClicks <- foldFun(AllClicksList) 
            folderName <- gsub('./', '', folder)
            write.csv(AllClicks, file=paste0('Taiki',folderName, '.CSV'))
        })
}
callFun <- function(x) {
    do.call(rbind, x)
}
microbenchmark(
    taikiDt = taikiTime(fastDt, rbindlist),
    taikiDf = taikiTime(fastDt, callFun),
    # jay = jayTime(),
    times = 1
)

microbenchmark(
    list_list = taikiTime(fastDt, rbindlist),
    list_call = taikiTime(fastDt, callFun),
    call_list = taikiTime(fastDf, rbindlist),
    call_call = taikiTime(fastDf, callFun),
    jayTime(),
    times=10
)

tf <- tempfile()
Rprof(tf)
taikiTime(fastDf)
# taikiTime(fastDt)
# jayTime()
Rprof(NULL)
summaryRprof(tf)

