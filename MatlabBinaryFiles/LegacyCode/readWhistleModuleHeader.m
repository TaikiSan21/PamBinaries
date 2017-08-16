function header = readModuleHeader(file)
header.length = fread(file, 1, 'int32');
header.identifier = fread(file, 1, 'int32');
header.version = fread(file, 1, 'int32');
header.binaryLength = fread(file, 1, 'int32');
% fseek(file, header.binaryLength, 'cof');
header.delayScale = fread(file, 1, 'int32');
fseek(file, 0, 'cof');