context('Test dataframe conversion')

test_that('PamBinaries objects convert to dataframe properly', {
    clickFile <- system.file('extdata', 'Click.pgdf', package='PamBinaries')
    clickData <- loadPamguardBinaryFile(clickFile)
    clickDf <- pbToDf(clickData)
    clickInfo <- loadPamguardBinaryFile(clickFile, skipData=TRUE)
    
    expect_equal(nrow(clickDf), 36)
    df <- expect_warning(pbToDf(clickData$data))
    expect_identical(clickDf, pbToDf(clickData))
    expect_identical(clickDf, df)
    expect_is(clickDf, 'data.frame')
    expect_null(pbToDf(clickInfo))
})

test_that('PamBinaries with match classifier templates', {
    clickFile <- system.file('extdata', 'ClickTemplate.pgdf', package='PamBinaries')
    clickData <- loadPamguardBinaryFile(clickFile)
    noTemplate <- pbToDf(clickData)
    clickDf <- pbToDf(clickData, templateNames = 1:6)
    
    expect_error(pbToDf(clickData, templateNames = 1:4))
    expect_warning(pbToDf(clickData, templateNames = 1:7))
    expect_is(clickDf, 'data.frame')
    expect_equal(ncol(clickDf), ncol(noTemplate) + 3 * 6)
    expect_identical(noTemplate, clickDf[, 1:ncol(noTemplate)])
})

test_that('PamBinaries with noise monitor data', {
    noiseFile <- system.file('extdata', 'NoiseMonitor.pgdf', package='PamBinaries')
    noiseData <- loadPamguardBinaryFile(noiseFile)
    noiseDf <- pbToDf(noiseData)
    
    expect_is(noiseDf, 'data.frame')
    expect_equal(nrow(noiseDf), 
                 length(noiseData$data) * noiseData$data[[1]]$nBands)
    expect_true(all(c('octaveBand', 'noiseMean', 'noisePeak') %in% colnames(noiseDf)))
})