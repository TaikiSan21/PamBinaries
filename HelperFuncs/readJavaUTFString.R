# Read a Java UTF-8 string. First 2 bytes are length of string, then the string.
# Based on MATLAB code by Michael Oswald

readJavaUTFString <- function(file) {
    # browser()
    len <- pamBinRead(file, 'int16', n=1)
    str <- pamBinRead(file, 'character', n=len)
    # print(len)
    list(len=len, str=str)
}