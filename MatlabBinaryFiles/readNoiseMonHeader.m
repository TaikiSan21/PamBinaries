function header = readNoiseMonHeader(file)
% reads module header information for the Noise Monitor module

header=readStdModuleHeader(file);
if (header.binaryLength~=0)
    header.nBands = fread(file, 1, 'int16');
    header.statsTypes = fread(file, 1, 'int16');
    header.loEdges = fread(file, header.nBands, 'float');
    header.hiEdges = fread(file, header.nBands, 'float');
end

