#' @title Read Pamguard Binary Data
#' 
#' @description A wrapper for reading various types of binary data.
#'   
#' @param fid The binary file being read
#' @param what The type of data to read. Int64 is not handled natively
#'   by R, see note.
#' @param n The number of objects to read.
#' @param seek Whether or not to just seek instead of reading
#' 
#' @return Data of the type and number specified.
#' 
#' @author Taiki Sakai \email{taiki.sakai@noaa.gov}
#' 
#' @note R does not natively support 64-bit integers. Current implementation
#'   is to read an int64 as 8 separate 1-byte raw pieces. These are converted
#'   from hexidecimal, shifted by the appropriate power of 2, then summed.
#'   Currently cannot read more than one int64 at a time, shouldn't be necessary.
#'   
pamBinRead <- function(fid, what=c('int8', 'int16', 'int32','int64',
                                   'float', 'character'), n, seek=FALSE) {
    endian <- 'big'
    if(seek) {
        switch(match.arg(what),
               int8 = seek(fid, n, origin='current'),
               int16 = seek(fid, n*2, origin='current'),
               int32 = seek(fid, n*4, origin='current'),
               int64 = seek(fid, n*8, origin='current'),
               float = seek(fid, n*4, origin='current'),
               character = seek(fid, n, origin='current'))
    } else if(n==0) {
        NA
    } else { switch(match.arg(what),
                    int8 = readBin(fid, 'integer', n=n, size=1, endian=endian),
                    int16 = readBin(fid, 'integer', n=n, size=2, endian=endian),
                    int32 = readBin(fid, 'integer', n=n, size=4, endian=endian),
                    int64 = {
                        if(n != 1) {
                            print(paste('WARNING: Currently does not accurately ',
                                        'read more than one int64 at a time.'))
                        }
                        # R has no int64, need to be creative. Reading in as raw - these are hexidecimal
                        a <- readBin(fid, 'raw', n=8, size=1, endian=endian)
                        # Convert each hex to int, then shift by power of 2
                        sum(strtoi(a, base=16)*2^(8*(7:0)))
                    },
                    float = readBin(fid, 'numeric', n=n, size=4, endian=endian),
                    character = rawToChar(readBin(fid, 'raw', n=n, size=1, endian=endian))
    )}
}
