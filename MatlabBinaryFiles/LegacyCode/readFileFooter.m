function footer = readFileFooter(file)

footer.length = fread(file, 1, 'int32');
footer.identifier = fread(file, 1, 'int32');
footer.nObjects = fread(file, 1, 'int32');
footer.dataDate = millisToDateNum(fread(file, 1, 'int64'));
footer.analysisDate = millisToDateNum(fread(file, 1, 'int64'));
footer.endSample = fread(file, 1, 'int64');
footer.fileLength = fread(file, 1, 'int64');
footer.endReason = fread(file, 1, 'int32');