# Pamguard binary reader helper
# Easier to translate from Java, uses bigendian
pamBinRead <- function(fid, what=c('int8', 'int16', 'int32',
                                   'int64', 'float', 'character'), n) {
    endian <- 'big'
    switch(match.arg(what),
           int8 = readBin(fid, 'integer', n=n, size=1, endian=endian),
           int16 = readBin(fid, 'integer', n=n, size=2, endian=endian),
           int32 = readBin(fid, 'integer', n=n, size=4, endian=endian),
           # int64 = {
           #     if(n != 1) {
           #         print(paste('WARNING: Currently does not accurately ',
           #                     'read more than one int64 at a time.'))
           #     }
           #     a <- as.numeric(readBin(fid, 'integer', n=n, size=4, endian=endian))
           #     b <- as.numeric(readBin(fid, 'integer', n=n, size=4, endian=endian))
           #     a*2^32 + b
           # },
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
           # character = readBin(fid, 'character', n=1, endian=endian)
    )
}
