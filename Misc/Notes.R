# NOTES #

# Dates can be converted with as.POSIXct(date, origin='1970-01-01', tz='UTC')

# TIMES: ClickTest: 11-12s vs 57s. 13955 events.
#        BigTest (DIFAR): 1.3s vs .5s 273 events.

# Possible future improvements - we should be able to read all of the data as one big chunk. 
# Either by looping through once just to get the data length, or can we work backwards from 
# the end of file? If we read in everything as a big raw chunk we can reorganize it into a
# matrix https://stat.ethz.ch/pipermail/r-help/2008-August/171645.html and then readBin from
# that structure all common types at a time. Would need to map the structure somehow. Possibly
# a matrix with type Id (stored as ints or factors) and length? Working backwards should work
# depending on what this 'datagram' nonsense is.

# Probably best to loop through once and store each data length and an event Id for it. Use this
# at end to split into separate events, then search through events for variable length data.

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
