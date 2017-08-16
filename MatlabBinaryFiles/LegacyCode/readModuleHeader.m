function header = readModuleHeader(file, readBinary)
if nargin < 2
    readBinary = false;
end
header.length = fread(file, 1, 'int32');
header.identifier = fread(file, 1, 'int32');
header.version = fread(file, 1, 'int32');
header.binaryLength = fread(file, 1, 'int32');
if readBinary
    header.binData = fread(file, header.binaryLength, 'int8');
else
    header.binData = [];
    fseek(file, header.binaryLength, 'cof');
end
fseek(file, 0, 'cof');