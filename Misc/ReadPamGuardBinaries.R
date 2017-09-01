######## ReadPamGuardBinaries.R     code from Taiki Sakai

# Package PamBinaries is on Github.  It is 25 separate R files so the easiest way to 
# get it is to install it with the devtools package.
# 
library(devtools)
devtools::install_github('TaikiSan21/PamBinaries')

library(PamBinaries)


# setwd("S:\\1650_PASCAL_2016\\Data\\DASBR_Pamguard_Post_Processing\\Binaries\\Final_Dataset")
setwd('./TestFiles/')
jayClicks <- function() {
    folders= list.dirs(recursive=F)
    nfolders= length(folders)
    ifolder=1     #use only for testing
    
    for (ifolder in 1:nfolders) {
        
        files= list.files(folders[ifolder],pattern=glob2rx("Click_Detector_Click_Detector_Clicks*.pgdf"),recursive=TRUE,full.names=TRUE)
        
        nfiles= length(files)
        for (ifile in 1:nfiles) {
            print(mem_used())
            cat(ifile," of ",nfiles,"  ",files[ifile],"\n")
            # Useage - loadPamguardBinaryFile will load any binary file type
            binaryData <- loadPamguardBinaryFile(files[ifile])
            
            # This is a list with two objects - data and fileInfo
            
            # To get your grouped click data:
            clickDf <- do.call(rbind, lapply(binaryData$data, function(x) {
                # Grabbing everything but the wave form
                waveId <- grep('wave', names(x))
                data.frame(x[-waveId])
            }))
            clickDf$ClickNo= 0:(length(clickDf$date)-1)
            clickDf$BinaryFile= sub(pattern="./",replacement="/",x=files[ifile])
            # This data is missing the columns in your output created by a 'clickAmplitude'
            # function that you had in Matlab. If you can send me that function I can add it.
            
            # To get waves together:
            nClicks <- length(binaryData$data)
            if (nClicks > 0) {
                waves <- vector('list', nClicks) # Initializing for speed
                for(i in 1:nClicks) {
                    waves[[i]] <- binaryData$data[[i]]$wave
                }
                if (ifile==1) {
                    AllClicks= clickDf
                } else {
                    AllClicks= rbind(AllClicks,clickDf) 
                }
            }
            
        }
        folderName= sub(pattern="./",replacement="",x=folders[ifolder])
        write.csv(AllClicks,file=paste(folderName,".CSV",sep=""))
        
    }
}
