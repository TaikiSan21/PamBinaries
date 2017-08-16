function header = readWMDHeader(file)
% reads module header information for the Whistle & Noise Detector module

header=readStdModuleHeader(file);
if (header.binaryLength~=0 && header.version>=1)
    header.delayScale = fread(file, 1, 'int32');
end

