# Clicks for J

#loadClickFile(fileName, arg2, arg3, arg4)

# arg2 is start date, arg3 is enddate, arg4 is preallocate.

matlabDate <- function(num) {
    as.POSIXct((num-719529)*86400, origin='1970-01-01', tz='UTC')
}
as.numeric(ymd_hms("2017-04-27 21:29:15 UTC")) -
as.numeric(ymd_hms("2017-03-09 04:26:28 UTC"))# me

as.numeric(ymd_hms("2014-08-06 18:41:47 UTC")) -
    as.numeric(ymd_hms("2014-06-18 01:39:00 UTC"))