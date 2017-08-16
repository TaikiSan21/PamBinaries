function header = readStdModuleHeader(file)
% reads the module header information common to all modules.  Differs from
% the legacy code in that it does not read in or skip any information
% specific to a module.  

header.length = fread(file, 1, 'int32');
header.identifier = fread(file, 1, 'int32');
header.version = fread(file, 1, 'int32');
header.binaryLength = fread(file, 1, 'int32');
