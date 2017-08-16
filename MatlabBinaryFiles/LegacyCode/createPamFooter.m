function footer = createPamFooter(nObjects, dataDate, endSample, fileLength )

footer.length = 8*4 + 4*4;
footer.identifier = -2;
footer.nObjects = nObjects;
footer.dataDate = dataDate;
footer.analysisDate = now();
footer.endSample = endSample;
footer.fileLength = fileLength;
footer.endReason = 0;