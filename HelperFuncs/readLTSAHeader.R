# Reads module header information for the LTSA module

readLTSAHeader <- function(file) {
      header <- readStdModuleHeader(file)
      if(header$binaryLength != 0) {
            header$fftLength <- pamBinRead(file, 'int32', n=1)
            header$fftHop <- pamBinRead(file, 'int32', n=1)
            header$intervalSeconds <- pamBinRead(file, 'int32', n=1)
      }
      return(header)
}