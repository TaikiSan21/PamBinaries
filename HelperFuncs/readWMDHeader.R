# Reads module header information for the Whistle & Moan Detector module

readWMDHeader <- function(file) {
    header <- readStdModuleHeader(file)
    if((header$binaryLength != 0) & (header$version >= 1)) {
        header$delayScale <- pamBinRead(file, 'int32', n=1)
    }
    return(header)
}