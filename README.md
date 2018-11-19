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

As of version 1.2 there is a second exported function, `pbToDf`. This function converts
the output of `loadPamguardBinaryFile` into a data frame. In order to do this it must
skip some of the information in the binaries that are not an appropriate size to put in
a data frame, as of version 1.2 the skipped data includes: any annotations data, click
wave forms, DIFAR demux data, and contour information from the Whistle & Moan Detector.
This function can be called directly, and is also called when as.data.frame is called 
on an object of class `PamBinary`; the output of `loadPamguardBinaryFile` is both `list`
and `PamBinary` as of version 1.2. 

The following should return identical results:

```r
df1 <- pbToDf(binaryData)
df2 <- as.data.frame(binaryData)
df3 <- data.frame(binaryData)
df4 <- pbToDf(binaryData$data)
```

The final example (`df4`) is not the preferred way to make this conversion, but is
included to avoid possibly surprising behavior. The function will issue a warning,
but will convert to a data frame as expected.

### Compatibility

PamBinaries should be compatible with Pamguard v2.00.15 and earlier.

#### Version 1.2.2

* Added support for multiple matched classifier annotations

#### Version 1.2

* The output of `loadPamguardBinaryFile` is now a `PamBinary` in addition to `list`.
This is just to allow for some easy method dispatch (like `as.data.frame` below),
and should not change any other behavior.

* Added a new exported function `pbToDf` that converts a `PamBinary` object to a
data frame. 

* Extended the generic `as.data.frame` for class `PamBinary`, `as.data.frame.PamBinary`
is just a wrapper for `pbToDf` but allows for quick and easy conversion to data frames.

#### Version 1.1

* Slight change to function output, conversion of date to POSIXct object is no longer
done automatically, only if the new flag `convertDate` is set to `TRUE`. Date is now
reported as seconds since 1970-01-01 UTC, and there is a second exported function
`convertPgDate` that will properly convert the numeric value to POSIXct. This change
was made for speed purposes, `loadPamguardBinaryFile` now runs approximately 40% faster
without the date conversion.

* Added support for Click Trigger Background binary files.

#### Version 1.0

* Initial release, see tutorial above.

### TO DO:

Add rest of annotation readers. Will probably need to re-adjust existing ones in future
PG versions.
