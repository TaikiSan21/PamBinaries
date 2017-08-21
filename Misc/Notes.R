# NOTES #

# Java UTF string wasn't behaving properly. Changed from int16 to int8 and it works
############## ^^^^^ ######### Stopped working with int8, now it works in int16
# Fuck you, R

# R has no 64 bit integers. Read in as two 32s, turned to numeric, then 
# a*2^32 + b. Seems to work fine, possibly sketchy.
# Fuck you, R

# Dates can be converted with as.POSIXct(date, origin='1970-01-01', tz='UTC')

# Dates arent working sometimes. The milliseconds are sometimes off by 2^32
# Probably because of my int64. Lets try to read in as raw. It reads one byte
# at a time as a hex. We can convert the hex to decimal, then shift those by
# 2^56, 2^48, ... 2^8, 2^0 or whatever power I need. strtoi(hex, base=16)

as.numeric(ymd_hms("2017-04-27 21:29:15 UTC")) -
    as.numeric(ymd_hms("2017-03-09 04:26:28 UTC")) #mine

as.numeric(ymd_hms("2014-08-06 18:41:47 UTC")) -
    as.numeric(ymd_hms("2014-06-18 01:39:00 UTC")) #mine
1407350507811.00 - a #matlab millis
a <- clicktest$dataSet[[1]]$millis
sprintf('%f', a)

c <- bigtest$dataSet[[1]]$millis/1000
sprintf('%f', c)


b <- test$dataSet[[1]]$millis/1000
sprintf('%f', b)
