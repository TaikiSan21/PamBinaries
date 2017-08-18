# NOTES #

# Java UTF string wasn't behaving properly. Changed from int16 to int8 and it works
############## ^^^^^ ######### Stopped working with int8, now it works in int16
# Fuck you, R

# R has no 64 bit integers. Read in as two 32s, turned to numeric, then 
# a*2^32 + b. Seems to work fine, possibly sketchy.
# Fuck you, R

# Dates can be converted with as.POSIXct(date, origin='1970-01-01', tz='UTC')