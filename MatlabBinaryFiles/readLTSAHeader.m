function header = readLTSAHeader(file)
% reads module header information for the LTSA module

header=readStdModuleHeader(file);
if (header.binaryLength~=0)
    header.fftLength = fread(file, 1, 'int32');
    header.fftHop = fread(file, 1, 'int32');
    header.intervalSeconds = fread(file, 1, 'int32');
end
