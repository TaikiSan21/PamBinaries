function footer = readFileFooterInfo(fid,version)
% reads in the binary file footer.  The input variable version is the file
% format read in from the file header.  As of version 3, the file footer
% includes the lowest and highest UID values in the file

footer.length = fread(fid, 1, 'int32');
footer.identifier = fread(fid, 1, 'int32');
footer.nObjects = fread(fid, 1, 'int32');
footer.dataDate = millisToDateNum(fread(fid, 1, 'int64'));
footer.analysisDate = millisToDateNum(fread(fid, 1, 'int64'));
footer.endSample = fread(fid, 1, 'int64');
if (version>=3)
    footer.lowestUID = fread(fid, 1, 'int64');
    footer.highestUID = fread(fid, 1, 'int64');
end
footer.fileLength = fread(fid, 1, 'int64');
footer.endReason = fread(fid, 1, 'int32');