# NOTES #

# Java UTF string wasn't behaving properly. Changed from int16 to int8 and it works
# R is somehow reading in the correct number of characters without specifying the length??
# Just having n=1 it reads the correct number. Seems sketchy. Can't find a way to have
# it read single char at a time like. Might be able to use raw(), looks like it reads in
# the chars as a code, could then convert them. 

# Header object length is not matching the length when stepping through.
#
# Need to compare to Matlab ftell going step by step.