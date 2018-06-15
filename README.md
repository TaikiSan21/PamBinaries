# PamBinaries

The PamBinaries package is a set of functions for reading PAMGuard binary output
(.pgdf files) into R. It is based on MATLAB functions created by Doug Gillespie
and Michael Oswald, and is written to work in the same way.

### Installation

Install the latest version from GitHub:

```r
# make sure you have Rtools installed
if(!require('devtools')) install.packages('devtools')
# install from GitHub
devtools::install_github('TaikiSan21/PamBinaries')

```

### Tutorial

The main function of this package is `loadPamguardBinaryFile`.

To read any PAMGuard binary file, just point this function to its location:

```r
myBinaryFile <- './Files/SomeBinaryFile.pgdf'
binaryData <- loadPamguardBinaryFile(myBinaryFile)
```
`loadPamguardBinaryFile` will automatically determine the type of data
stored in the binary file (Click Detector, DIFAR Processing, etc.) and call
the appropriate function. 

`binaryData` is a list with two objects, `data`
and `fileInfo`. `binaryData$data` is a list with one element for every 
individual item in the binary file, the contents of which will be different
depending on the type of binary file being read. `binaryData$fileInfo` 
is a list containing some metadata and information used by the 
`loadPamguardBinaryFile` function.

There is also an option to reduce the size of objects read from a 
Click Detector, Whistle and Moan Detector, or DIFAR:

```r
binaryDataSmall <- loadPamguardBinaryFile(myBinaryFile, skipLarge=TRUE)
```

This option will skip reading the waveform for click and difar data, and
skip the contour and slice information for whistle data. This greatly reduces
the size of the loaded object and the time required to read the data, useful
if you need to process a large number of files.

You can also choose to only read in data with specific UIDs:

```r
binarySomeUIDs <- loadPamguardBinaryFile(myBinaryFile, keepUIDs=c(100000001, 100000002))
```

This only read in the data with the UIDs provided in the keepUIDs argument, skipping over
the rest. This can greatly reduce the time required to read data if you only need a
small subset of the items in the binary file.

### TO DO:




