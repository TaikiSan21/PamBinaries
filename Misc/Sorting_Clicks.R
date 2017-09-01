dir.create("C:/Griffiths/PASCAL/SpeciesSpecific/NBHF/binaryextraction_output")
setwd("C:/Griffiths/PASCAL/SpeciesSpecific/NBHF/binaryextraction_output")


#Get all your packages
library(seewave)
library(RSQLite)
library(plyr)
library(tuneR)
library(signal)

#Read in the database
sqlite <- dbDriver("SQLite")
StaDB <- dbConnect(sqlite,"F:/PASCAL/Data/PamGuard/Database/NBHF/NBHF Station-2 ST4300-B Ch1trigger.sqlite3")
allevents <- dbReadTable(StaDB, "Click_Detector_OfflineEvents")
allclicks <- dbReadTable(StaDB, "Click_Detector_OfflineClicks")

#Read in the Binary output
bincsv=list.files("H:/Sta2/", pattern= "click", recursive=TRUE,full.names = TRUE)
wavecsv=list.files("H:/Sta2/", pattern= "wavef",recursive=TRUE, full.names = TRUE)

bin=NULL
for (i in seq_along(bincsv)) {
  a=read.csv(bincsv[i], stringsAsFactors = FALSE)
  a$binDT=substr(bincsv[i],25,39)
  a$datetime=as.POSIXct(a$millis/1000, origin = "1970-01-01", tz="UTC")
  bin=rbind(bin,a)
}

#Recreate click number to start at 0 like PG likes it.
bin$ClickNo <- (ave(bin$angleErrors,bin$binDT, FUN = seq_along))-1

wave=NULL
for (i in seq_along(wavecsv)) {
  a=read.csv(wavecsv[i], header = FALSE)
  a$binDT=substr(wavecsv[i],25,39)
  wave=rbind(wave,a)
}


#Divide the wave into individual clicks, and other imporant varabables.
binsize=5000
NumChs=2
f=288000
FFTsize=288
station="Sta2"

#Split by binary time stamp
wavebins=split.data.frame(wave,wave$binDT)

#From the database, list all of binary files that have events present.
events=unique(allclicks$BinaryFile)
eventDT=substr(events,38,52)
#For later, files and corresponding event numbers
eventID=as.data.frame(cbind(allclicks$EventId,allclicks$BinaryFile))
eventID=eventID[!duplicated(eventID),]

#Subset the data based on the binary list
x=unlist(names(wavebins))
i1 <- x %in% eventDT
waveevents=wavebins[i1]


#split each event into clicks by bin size
waveforms=list()
for(i in seq_along(waveevents)) {
  tempobj=split(waveevents[[i]]$V1,ceiling(seq_along(waveevents[[i]]$V1)/(binsize*NumChs)))
  name= names(waveevents)[[i]]
  waveforms[[name]]=tempobj
}


#split each click into ch0 and ch1.
ch0=list()
for (i in seq_along(waveforms)) {
  a=list()#names(waveforms)[[i]])
  b=for (j in seq_along(waveforms[[i]])) {
    tempobj=head((waveforms[[i]])[[j]],n=binsize)
    name <- paste('click',seq_along((waveforms)[[i]])[[j]]-1,sep='')
    a[[name]] <- tempobj
  }
  name1 <- names(waveforms)[[i]]
  ch0[[name1]] <- a
}

ch1=list()
for (i in seq_along(waveforms)) {
  a=list()
  b=for (j in seq_along(waveforms[[i]])) {
    tempobj=tail((waveforms[[i]])[[j]],n=binsize)
    name <- paste('click',seq_along((waveforms)[[i]])[[j]]-1,sep='')
    a[[name]] <- tempobj
  }
  name1 <- names(waveforms)[[i]]
  ch1[[name1]] <- a
}

#include clicks only selected for events. 

#Create a UID for the clicks in the database.
allclicks$click=paste0("Click",allclicks$ClickNo)
allclicks$EventId=paste0("Event",allclicks$EventId)
allclicks$binDT=substr(allclicks$BinaryFile, 38,52)
allclicks$UiD=paste0(station,"_",allclicks$EventId,"_",allclicks$binDT,"_",allclicks$click)
allclicks$sUiD=paste0(station,"_",allclicks$binDT,"_",allclicks$click)

#create a temporary UID for all of the clicks in the binaries.
bin$click=paste0("Click",bin$ClickNo)
bin$UiD=paste0(station,"_UNKN_",bin$binDT,"_",bin$click)
bin$sUiD=paste0(station,"_",bin$binDT,"_",bin$click)

#match and filter based on the UiD.  Get a list of clicks excluded.
#data=NULL
#for(d in 1:nrow(allclicks)) {
#  k<-unique(substr(allclicks$UiD,0,4)) == unique(substr(bin$UiD,0,4))
#  if(k==FALSE) stop("Project ID is not identical.  Make sure your data is from the same PG database")
#}

data=merge(allclicks,bin, by="sUiD", all = TRUE)
#Create a UiD column where clicks that have events are a full UiD code, and clicks not in events are marked as "UNKN"
data$UiD=ifelse(is.na(data$UiD.x),data$UiD.y,data$UiD.x)
justevts=data[!is.na(data$EventId),]

#Like the waveform data, split it by the Datetime stamp
byDT=split.data.frame(data,data$binDT.y)

#Generate a list that is strcutred like the waveform data
sortUNKN=list()
for (i in seq_along(byDT)) {
  a=as.list(byDT[[i]]$UiD)
  names(a)=byDT[[i]]$click.y
  n=names(byDT[i])
  sortUNKN[[n]]<-a
}

#Mark the clicks that you want to include in further analysis (clicks in events)
rmclicks=lapply(sortUNKN, function(x) lapply(x, function(y) grep("Event",y)))
i1 <- lapply(rmclicks,function(x) sapply(x,any))

#Subset the clicks from the waveform list.
SelectClicks=list()
for(j in seq_along(waveforms)) {
  a=waveforms[[j]][i1[[j]]]
  name <- names(waveforms)[[j]]
  SelectClicks[[name]] <- a
}

#Split into channels
ch0=list()
for (i in seq_along(SelectClicks)) {
  a=list()
  b=for (j in seq_along(SelectClicks[[i]])) {
    tempobj=head((SelectClicks[[i]])[[j]],n=binsize)
    name <- names(SelectClicks[[i]])[[j]]
    a[[name]] <- tempobj
  }
  name1 <- names(SelectClicks)[[i]]
  ch0[[name1]] <- a
}

ch1=list()
for (i in seq_along(SelectClicks)) {
  a=list()
  b=for (j in seq_along(SelectClicks[[i]])) {
    tempobj=tail((SelectClicks[[i]])[[j]],n=binsize)
    name <- names(SelectClicks[[i]])[[j]]
    a[[name]] <- tempobj
  }
  name1 <- names(SelectClicks)[[i]]
  ch1[[name1]] <- a
}

#Save your lists
dput(ch0, file = "Station2_Channel0.txt")
dput(ch1, file = "Station2_Channel1.txt") 

#To retrieve: dget()

name=names(justevts)
cleandata=justevts[,c(1:26,28)]
name=name[c(1:26,28)]
name=gsub(".x","",name)
names(cleandata)=name

write.csv(cleandata, "Click_binaryNoffline_info.csv", row.names = FALSE)
