# Count the number of set bits in the channel map
# Changed from 32 -> 30 because R hits the maximum integer limit
# There shouldn't ever be that many channels.
countChannels <- function(channelMap) {
    nC <- 0
    j <- 1
    for(i in 1:30) { # Changed, see above
        if(bitwAnd(channelMap, j) != 0) {
            nC <- nC + 1
        }
        j <- j * 2
    }
    return(nC)
}
