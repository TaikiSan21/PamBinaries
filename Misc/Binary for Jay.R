# Package is on Github.  It is 25 separate R files so the easiest way to 
# get it is to install it with the devtools package.
# Let me know if this doesn't work.

# I think Shannon wanted to check some things when she got back before we 
# start sharing this with everyone.
library(devtools)
devtools::install_github('TaikiSan21/PamBinaries')

library(PamBinaries)

# Useage - loadPamguardBinaryFile will load any binary file type
filePath <- './TestFiles/Click_Detector_Click_Detector_Clicks_20160903_004214.pgdf'
binaryData <- loadPamguardBinaryFile(filePath)

# This is a list with two objects - data and fileInfo

# To get your grouped click data:
clickDf <- do.call(rbind, lapply(binaryData$data, function(x) {
    # Grabbing everything but the wave form
    waveId <- grep('wave', names(x))
    data.frame(x[-waveId])
}))
# This data is missing the columns in your output created by a 'clickAmplitude'
# function that you had in Matlab. If you can send me that function I can add it.

# To get waves together:
nClicks <- length(binaryData$data)
waves <- vector('list', nClicks) # Initializing for speed
for(i in 1:nClicks) {
    waves[[i]] <- binaryData$data[[i]]$wave
}
tmp <- tempfile()
Rprof(tmp)
# binaryData <- loadPamguardBinaryFile(filePath)
clickDf <- do.call(rbind, lapply(binaryData$data, function(x) {
    # Grabbing everything but the wave form
    waveId <- grep('wave', names(x))
    data.frame(x[-waveId])
}))
Rprof(NULL)

summaryRprof(tmp)


tmp <- tempfile()
Rprof(tmp)
# binaryData <- loadPamguardBinaryFile(filePath)
clickDf2 <- new(binaryData$data)
Rprof(NULL)

summaryRprof(tmp)


library(microbenchmark)
reg <- function(data) {
    do.call(rbind, lapply(data, function(x) {
        waveId <- grep('wave', names(x))
        data.frame(x[-waveId])
    }))}

new <- function(data) {
    do.call(rbind, lapply(data, function(x) {
        tooBig <- which(sapply(x, function(y) !is.null(dim(y))))
        x <- x[-tooBig]
        class(x) <- 'data.frame'
        attr(x, 'row.names') <- .set_row_names(length(x[[1]]))
        x
    }))}

microbenchmark(
    reg = reg(um),
    new = new(um),
    times=10
)

t <- binaryData$data[[1]]
t <- t[-waveId]
class(t) <- 'data.frame'
attr(t, "row.names") <- .set_row_names(length(t[[1]]))

quickdf <- function(l) {
    # skip <- grep('b', names(l))
    skip <- which(sapply(l, function(y) !is.null(dim(y))))
    skip <- 24
    l <- l[-skip]
    class(l) <- "data.frame"
    attr(l, "row.names") <- .set_row_names(length(l[[1]]))
    l
}
slowdf <- function(l) {
    # skip <- grep('b', names(l))
    skip <- which(sapply(l, function(y) !is.null(dim(y))))
    skip <- 24
    data.frame(l[-skip])
}

l <- lapply(1:26, function(i) runif(1e3))
names(l) <- letters

microbenchmark(
    quick_df      = quickdf(l),
    slow_df = slowdf(l),
    as.data.frame = data.frame(l[-24]),
    times=1000
)
