# Pamguard binary reader helper
# Easier to translate from Java, uses bigendian
pamBinRead <- function(fid, what=c( 'int8', 'int16', 'int32', 'int64', 'float', 'character'), n) {
      endian <- 'big'
      switch(match.arg(what),
             int8 = readBin(fid, 'integer', n=n, size=1, endian=endian),
             int16 = readBin(fid, 'integer', n=n, size=2, endian=endian),
             int32 = readBin(fid, 'integer', n=n, size=4, endian=endian),
             int64 = readBin(fid, 'integer', n=n, size=8, endian=endian),
             float = readBin(fid, 'numeric', n=n, size=4, endian=endian),
             character = readBin(fid, 'character', n=n, endian=endian))
}
      