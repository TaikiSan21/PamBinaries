context('test Pamguard binary reading and basic unit conversion')

test_that('Pamguard binary files are able to be read', {
    clickFile <- system.file('extdata', 'Click.pgdf', package='PamBinaries')
    clickData <- loadPamguardBinaryFile(clickFile)
    clickInfo <- loadPamguardBinaryFile(clickFile, skipData = TRUE)
    
    expect_s3_class(clickData, 'PamBinary')
    expect_identical(clickData$fileInfo, clickInfo$fileInfo)
    expect_equal(length(clickData$data), 36)
    expect_equal(length(clickInfo$data), 0)
    noFile <- expect_warning(loadPamguardBinaryFile('FILEDNE'))
    expect_null(noFile)
})

test_that('Converting units of PamBinary data', {
    wmFile <- system.file('extdata', 'WM.pgdf', package='PamBinaries')
    wmData <- loadPamguardBinaryFile(wmFile)
    wmInfo <- loadPamguardBinaryFile(wmFile, skipData=TRUE)
    wmFreq <- contourToFreq(wmData)
    clickFile <- system.file('extdata', 'Click.pgdf', package='PamBinaries')
    clickData <- loadPamguardBinaryFile(clickFile)
    expect_equal(wmFreq$data[[1]]$freq[1],
                 wmFreq$data[[1]]$contour[1] * 44100 / 2048)
    expect_true(all(c('freq', 'allFreq', 'time') %in% names(wmFreq$data[[1]])))
    expect_error(contourToFreq(clickData))
    expect_equal(length(contourToFreq(wmInfo)$data), 0)
    
    expect_equal(convertPgDate(1559394732),
                 as.POSIXct("2019-06-01 13:12:12 UTC", tz='UTC'))
    expect_true(is.na(convertPgDate(NA)))
    expect_null(convertPgDate(NULL))
})
    