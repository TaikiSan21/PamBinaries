# PamBinaries 1.8.1

* `contourToFreq` sometimes not working for all files, better method now

# PamBinaries 1.8.0

* New Gemini Track Reader code added from pull request

# PamBinaries 1.7.1

* `keepUIDs` option logic changed to be more consistent in the face of strangeness 

# PamBinaries 1.7.0

* Adding support for deep learning module

* Added some missing annotation functions

# PamBinaries 1.6.3

* Odd issues with `loadBackgroundNoise` not passing tests

# PamBinaries 1.6.2

* `pbToDf` was failing on GPL data

# PamBinaries 1.6.1

* `pbToDf` was failing on some binary types with a `$noise` value that were
not the Noise Band Monitor

# PamBinaries 1.6.0

* Adding functions `loadBackgroundNoise`, `plotBackgroundNoise`, and
`combineBackgroundNoise` to assist with
working with new background noise data in binary files

* Added new test/example files for background noise

# PamBinaries 1.5.1

* Fixed issue when reading GPL with empty contours

# PamBinaries 1.5.0

* Michael Oswald added support for GPL detector

# PamBinaries 1.4.1

* Fixed issue with noise band monitor data in `pbToDf` related to `vctrs` package update

# PamBinaries 1.4.0

* Added a `NEWS.md` file to track changes to the package.

* Added examples to all functions in prep for CRAN submission

* Added unit testing

* Added files to inst/extdata for examples and testing

* Fixed a few functions to fail more gracefully on bad/NA/NULL input

* `pbToDf` either fails or warns on `templateNames` number mismatch depnding
on whether too few or too many names were provided

* Added support for Travis-CI

# PamBinaries 1.3.5

* `pbToDf` now works for noise band monitor data in a well-behaved manner

# PamBinaries 1.3.4

* `pbToDf` now reads match and reject corr for click template classifier, and also fixed
an issue where it would crash sadly if templateNames were supplied but a click did not have
a corresponding number of templates

# PamBinaries 1.3.3

* Clicks now also store `maxAmplitude` in the output

# PamBinaries 1.3.2

* Added feature to `pbToDf` so read click template classifier thresholds

# PamBinaries 1.3.1

* Fixed `pbToDf` so that it is much faster (~50x)

# PamBinaries 1.3.0

* Added option `skipData` to `loadPamguardBinaryFile` that will only read in file headers
and footers. This also speeds up data loads when using `keepUIDs` argument, and preserves
file footers in this case (did not before)

* `skipLarge` was actually slower with Whistle and Moan data, now is much faster and does
not save any of the contour data to the data file (it was all 0s before, so not meaningful)

# PamBinaries 1.2.7

* Better `plotWMD` labeling and `verbose` option for `contourToFreq` that prints parameters

# PamBinaries 1.2.6

* Fixed a bug in SR and FFT parameter calculation for `contourToFreq`

# PamBinaries 1.2.5

* Added function `plotWMD` to look at whistle contour plots

# PamBinaries 1.2.4 

* Added function `contourToFreq` for adding frequency and time information to WMD binaries

# PamBinaries 1.2.3

* Fixed typo in matched classifier annotations

# PamBinaries 1.2.2

* Added support for multiple matched classifier annotations

# PamBinaries 1.2

* The output of `loadPamguardBinaryFile` is now a `PamBinary` in addition to `list`.
This is just to allow for some easy method dispatch (like `as.data.frame` below),
and should not change any other behavior.

* Added a new exported function `pbToDf` that converts a `PamBinary` object to a
data frame. 

* Extended the generic `as.data.frame` for class `PamBinary`, `as.data.frame.PamBinary`
is just a wrapper for `pbToDf` but allows for quick and easy conversion to data frames.

# PamBinaries 1.1

* Slight change to function output, conversion of date to POSIXct object is no longer
done automatically, only if the new flag `convertDate` is set to `TRUE`. Date is now
reported as seconds since 1970-01-01 UTC, and there is a second exported function
`convertPgDate` that will properly convert the numeric value to POSIXct. This change
was made for speed purposes, `loadPamguardBinaryFile` now runs approximately 40% faster
without the date conversion.

* Added support for Click Trigger Background binary files.

# PamBinaries 1.0

* Initial release, see tutorial above.
