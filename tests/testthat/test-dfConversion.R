context('Test dataframe conversion')

test_that('PamBinaries objects convert to dataframe properly', {
    clickFile <- system.file('extdata', 'Click.pgdf', package='PamBinaries')
    clickData <- loadPamguardBinaryFile(clickFile)
    clickDf <- pbToDf(clickData)
    clickInfo <- loadPamguardBinaryFile(clickFile, skipData=TRUE)
    
    expect_equal(nrow(clickDf), 36)
    expect_identical(clickDf, pbToDf(clickData$data))
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