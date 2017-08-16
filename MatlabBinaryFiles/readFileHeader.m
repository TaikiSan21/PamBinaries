function header = readFileHeader(file, readExtra)
if nargin < 2
    readExtra = 0;
end
header.length = fread(file, 1, 'int32');
header.identifier = fread(file, 1, 'int32');
header.fileFormat = fread(file, 1, 'int32');
header.pamguard = char(fread(file, 12, 'uchar')');
header.version = readJavaUTFString(file);
header.branch = readJavaUTFString(file);
header.dataDate = millisToDateNum(fread(file, 1, 'int64'));
header.analysisDate = millisToDateNum(fread(file, 1, 'int64'));
header.startSample = fread(file, 1, 'int64');
header.moduleType = readJavaUTFString(file);
header.moduleName = readJavaUTFString(file);
header.streamName = readJavaUTFString(file);
header.extraInfoLen = fread(file, 1, 'int32');
if readExtra
    header.extraInfo = fread(file, header.extraInfoLen, 'int8');
else
    fseek(file, header.extraInfoLen, 'cof');
end
