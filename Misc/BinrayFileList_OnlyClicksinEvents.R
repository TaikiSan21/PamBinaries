setwd("C:/Users/emily.griffiths/Documents/PASCAL/")


# You'll need to install my PamBinaries package (I updated it recently)
# library(devtools)
# devtools::install_github('TaikiSan21/PamBinaries')

# Change this path to wherever you want to put this file
source('clickCombine.R')

library(dplyr)
library(RSQLite)

#load BW databases.  These are the databases for BWs only, all other species events have been deleted.
filelist=list.files("S:/1650_PASCAL_2016/Data/DASBR_Pamguard_Post_Processing/Database/Final_BWDataset_20170601", pattern = "sqlite3", full.names = TRUE, recursive = TRUE)

sqlite <- dbDriver("SQLite")
ClicksNEvents=NULL
binlist=NULL
for (i in 1:length(filelist)) {
  conn <- dbConnect(sqlite,filelist[i])
  eventtable <- dbReadTable(conn, "Click_Detector_OfflineEvents")
  clicktable <- dbReadTable(conn, "Click_Detector_OfflineClicks")
  CNE=merge(clicktable,eventtable, by.x="EventId", by.y="Id",all.y = TRUE)
  shorthand=CNE[,c("EventId","BinaryFile","ClickNo")]
  temp=basename(filelist[i])
  station=gsub( "_MAS.*$", "", temp)
  shorthand$Station=station
  path=paste0("S:/1650_PASCAL_2016/Data/DASBR_Pamguard_Post_Processing/Binaries/Full_Click_Runs/Processed_Master_Copies/Binaries_",station)
  binlist_short=cbind(shorthand["BinaryFile"],path, station)
  ClicksNEvents=rbind(ClicksNEvents,CNE)
  binlist=rbind(binlist,binlist_short)
}
###### LOADING WAVE FILES FROM BINARIES ########

# Need to clean up blank space for use later
ClicksNEvents$BinaryFile <- gsub(' ', '', ClicksNEvents$BinaryFile) 
binlist$BinaryFile <- gsub(' ', '', binlist$BinaryFile)

# Pick whatever event you want to look at, then filter down the main dataframe
myEvent <- 23
myClicksNEvents <- filter(ClicksNEvents, EventId==myEvent)

# This is where we'll store the waves
waves <- list()

# For each unique binary file in our event...
for(file in unique(myClicksNEvents$BinaryFile)) {
    # There were some NA entries so I put this here for safety. Might want to look into why they are NA
    if(is.na(file)) {
        print('WARNING: Binary file iS NA')
        next
    }
    ####### THIS PART MIGHT BREAK #######
    # I need to get the exact path to each binary file, which is something like
    # .../Station/20160807/BinaryFile.pgdf
    # The date numbers in the middle I'm getting from part of the name of the binary file
    # Click_Detector_....._20160807_1234.pgdf
    # If that isn't how that folder is always named, this part will break.
    subFolder <- gsub('.*_([0-9]*)_[0-9]*\\.pgdf$', '\\1', file)
    ######################################
    
    # Full file path to the binary we want
    filePath <- paste(filter(binlist, BinaryFile==file)$path[1], subFolder, file, sep='/')
    # Get a list of every ClickNo in this event for this binary file
    clickNos <- filter(myClicksNEvents, BinaryFile==file)$ClickNo
    # Read binary
    waves[file] <- combineClickFiles(fileList=filePath, getWave=TRUE, 
                                                clickNos = clickNos)
}
# waves is a list with 1 element for each unique binary file in that event, and each of these
# is named whatever the file name is. Each file has a list of every click you want. So if you
# need to get ClickNo 52 from file BinaryFile1.pgdf you can just do
# waves$BinaryFile1.pgdf$`52`
    


write.csv(ClicksNEvents, "BW Clicks in Events PASCAL.csv", row.names = FALSE)
write.csv(ClicksNEvents, "BW Binary File List.csv", row.names = FALSE)
