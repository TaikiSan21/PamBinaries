# Comparison to Jay's
library(data.table)
library(dplyr)
library(microbenchmark)
library(R.utils)
library(pryr)
sourceDirectory('./R')
source('./Misc/clickCombine.R')

microbenchmark(
    taiki=taikiClicks(quiet=TRUE, getWave=FALSE),
    # allcheck=taikiClicks2(quiet=TRUE),
    jay = jayClicks(),
    times=1)

try <- 'S:\\1650_PASCAL_2016\\Data\\DASBR_Pamguard_Post_Processing\\Binaries\\Full_Click_Runs\\Processed_Master_Copies\\Binaries_Station-23_Soundtrap-E'
taikiClicks(folders=try, getWave=FALSE)

trySmall <- 'S:\\1650_PASCAL_2016\\Data\\DASBR_Pamguard_Post_Processing\\Binaries\\Full_Click_Runs\\Processed_Master_Copies\\Binaries_Station-8_Card-B'
taikiClicks(folders=trySmall, getWave=FALSE)


taikiClicks(folders='./TestFiles/1', getWave=FALSE)
test <- loadPamguardBinaryFile('./TestFiles/ClickTest.pgdf')

microbenchmark(
    wave = loadPamguardBinaryFile('./TestFiles/ClickTest.pgdf', getWave=TRUE),
    noWave = loadPamguardBinaryFile('./TestFiles/ClickTest.pgdf', getWave=FALSE),
    times=5
)

microbenchmark(
    wave = taikiClicks(dir='./TestFiles/', quiet=TRUE, getWave=TRUE),
    noWave = taikiClicks(dir='./TestFiles/', quiet=TRUE, getWave=FALSE),
    times=3
)

waveDf2 <- rbindlist(lapply(noWave$data, function(d) {
    class(d) <- 'data.frame'
    attr(d, 'row.names') <- .set_row_names(length(d[[1]]))
    d}))
