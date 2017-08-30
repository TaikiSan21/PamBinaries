# Difar testing

difar <- loadPamguardBinaryFile('./TestFiles/BigTest.pgdf')
demux <- difar$data[[1]]$demuxData

bscan <- difarBscan(demux[,1], demux[,2], demux[,3], 
                    fftLength=512, sampleRate=3000, freqRange=c(1000,1250), 
                    degstep=2, outputType = 'MVDR')

plot_ly(x= bscan$ang, y=bscan$freqs, z=bscan$Output, type='heatmap')