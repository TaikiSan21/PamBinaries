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
#' @author Taiki Sakai \email{taiki.sakai@@noaa.gov}
#' 
#' @note R does not natively support 64-bit integers. Current implementation
#'   is to read an int64 as 8 separate 1-byte raw pieces. These are converted
#'   from hexidecimal, shifted by the appropriate power of 2, then summed.
#'   Currently cannot read more than one int64 at a time, shouldn't be necessary.
#'   
pamBinRead <- function(fid, what=c('int8', 'int16', 'int32','int64', 'uint8', 'uint16',
                                   'float', 'double', 'character'), n, seek=FALSE) {
    endian <- 'big'
    
    if(length(n) == 0) {
        stop('"n" is length 0, it appears there is no data to read.')
    }
    
    if(seek) {
        return(
            switch(what,
                   int8 = seek(fid, n, origin='current'),
                   uint8 = seek(fid, n, origin='current'),
                   int16 = seek(fid, n*2, origin='current'),
                   uint16 = seek(fid, n*2, origin='current'),
                   int32 = seek(fid, n*4, origin='current'),
                   int64 = seek(fid, n*8, origin='current'),
                   float = seek(fid, n*4, origin='current'),
                   double = seek(fid, n*8, origin='current'),
                   character = seek(fid, n, origin='current'),
                   warning(paste0("Can't read binary data type ", what))
            )
        )
    }
    
    if(n == 0) {
        return(NA)
    }
    
    switch(what,
           int8 = readBin(fid, 'integer', n=n, size=1, endian=endian),
           uint8 = readBin(fid, 'integer', n=n, size=1, signed = FALSE, endian=endian),
           int16 = readBin(fid, 'integer', n=n, size=2, endian=endian),
           uint16 = readBin(fid, 'integer', n=n, size=2, signed = FALSE, endian=endian),
           int32 = readBin(fid, 'integer', n=n, size=4, endian=endian),
           int64 = {
               if(n != 1) {
                   print(paste('WARNING: Currently does not accurately ',
                               'read more than one int64 at a time.'))
               }
               # R has no int64, need to be creative. Reading in as raw - these are hexidecimal
               a <- readBin(fid, 'raw', n=8, size=1, endian=endian)
               # Convert each hex to int, then shift by power of 2
               a <- sum(strtoi(a, base=16)*2^(8*(7:0)))
               if(a <= .Machine$integer.max) {
                   as.integer(a)
               } else {
                   a
               }
           },
           float = readBin(fid, 'numeric', n=n, size=4, endian=endian),
           double = readBin(fid, 'numeric', n=n, size=8, endian=endian),
           character = rawToChar(readBin(fid, 'raw', n=n, size=1, endian=endian)),
           warning(paste0("Can't read binary data type ", what))
    )
}
