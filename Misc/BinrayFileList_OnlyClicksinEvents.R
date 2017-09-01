setwd("C:/Users/emily.griffiths/Documents/PASCAL/")


library(RSQLite)

#load BW databases.  These are the databases for BWs only, all other species events have been deleted.
filelist=list.files("S:/1650_PASCAL_2016/Data/DASBR_Pamguard_Post_Processing/Database/Final_Dataset", pattern = "sqlite3", full.names = TRUE, recursive = TRUE)

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

write.csv(ClicksNEvents, "BW Clicks in Events PASCAL.csv", row.names = FALSE)
write.csv(ClicksNEvents, "BW Binary File List.csv", row.names = FALSE)
