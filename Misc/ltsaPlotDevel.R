#ltsa plotting

# loads LTSA data from binary files into a plottable matrix
# INPUTS: 
#   ltsaFiles - binary folder of ltsa binaries, or a vector of the specific ones you want
#   vp2p - peak 2 peak voltage
#   sens - sensitivity of equipment
#   gain - gain of equipment
#   dateRange - range of times you want to plot, in POSIXct format c(start, end)
#   freqRange - frequency range you want to plot, in Hertz c(min, max)
# OUTPUTS:
#   A list with 3 components:
#    $ltsaMatrix - a matrix of the LTSA values (each column is a separate time, row is frequency)
#    $date - the dates corresponding to the columns
#    $freq - the frequencies (Hz) corresponding to the rows
loadLtsaMatrix <- function(ltsaFiles, vp2p=2, sens=20, gain=-201, dateRange=NULL, freqRange=NULL) {
    if(length(ltsaFiles) == 1 &&
       dir.exists(ltsaFiles)) {
        ltsaFiles <- list.files(ltsaFiles, '.*LTSA.*pgdf$')
    }
    ltsaOut <- lapply(ltsaFiles, function(x) {
        ltsa <- loadPamguardBinaryFile(x, convertDate=FALSE)
        
        fftLen <- ltsa$fileInfo$moduleHeader$fftLength
        sr <- ltsa$fileInfo$moduleHeader$fftHop / ltsa$fileInfo$moduleHeader$intervalSeconds
        binWidth <- sr / fftLen
        freq <- (1:(fftLen/2)) * sr / fftLen
        date <- sapply(ltsa$data, function(x) x$date)
        if(is.null(dateRange)) {
            inDate <- rep(TRUE, length(date))
        } else {
            inDate <- date >= dateRange[1] & date <= dateRange[2]
        }
        if(is.null(freqRange)) {
            inFreq <- rep(TRUE, length(freq))
        } else {
            inFreq <- freq >= freqRange[1] & freq <= freqRange[2]
        }
        ltsaMatrix <- matrix(0, 
                         nrow=sum(inFreq),
                         ncol=length(ltsa$data))
        for(i in seq_along(ltsa$data)) {
            ltsaMatrix[, i] <- ltsa$data[[i]]$data[inFreq]
        }
        ltsaMatrix <- ltsaMatrix[, inDate]
        ltsaMatrix <- ((ltsaMatrix / fftLen) * sqrt(2)) / sqrt(binWidth)
        ltsaMatrix <- 20*log10((ltsaMatrix / 2) * vp2p) - (sens + gain)
        list(date=date[inDate], freq=freq[inFreq], ltsaMatrix=ltsaMatrix)
    })
    date <- unlist(lapply(ltsaOut, function(x) x$date))
    date <- as.POSIXct(date, origin='1970-01-01', tz='UTC')
    freq <- lapply(ltsaOut, function(x) x$freq)
    if(length(freq) > 1 && !purrr::reduce(freq, .f=identical)) {
        stop('Not all frequency ranges are identical, cannot combine LTSAs')
    }
    freq <- freq[[1]]
    # browser()
    ltsaMatrix <- purrr::reduce(lapply(ltsaOut, function(x) x$ltsaMatrix), cbind)
    list(date=date, freq=freq, ltsaMatrix=ltsaMatrix)
}

library(PamBinaries)
ltsaFile <- '~/../Downloads/LTSA_Long_Term_Spectral_Average_LTSA_20191013_180005.pgdf'
ltsa <- loadLtsaMatrix(ltsaFile)
# matlab-like colors
jetColors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
image(x=as.numeric(ltsa$date), y=ltsa$freq, z=t(ltsa$ltsaMatrix), useRaster = TRUE, xaxt='n', col=jetColors(32),
      xlab='Date', ylab='Frequency (Hz)')
axis(side=1, labels=pretty(ltsa$date), at=pretty(as.numeric(ltsa$date)))
