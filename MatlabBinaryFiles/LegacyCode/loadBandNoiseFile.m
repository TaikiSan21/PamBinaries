function [noises fileHeader fileFooter moduleHeader moduleFooter] = loadBandNoiseFile(fileName, arg2, arg3, arg4)
noises = [];
fileHeader = [];
fileFooter = [];
moduleHeader = [];
moduleFooter = [];
if nargin == 0
    fileName = 'C:\WhistleClassifier\48kBS\19900101\WhistlesMoans_Whistle_and_Moan_Detector_Contours_19900101_001000.pgdf'
end
if nargin < 2
    preallocate = 0;
    startDate = 0;
    endDate = datenum('1 Jan 2099');
end
if nargin == 2
    preallocate = arg2;
    startDate = 0;
    endDate = datenum('1 Jan 2099');
end
if nargin == 3
    preallocate = 0;
    startDate = arg2;
    endDate = arg3;
end
if nargin == 4
    preallocate = arg4;
    startDate = arg2;
    endDate = arg3;
end
try
f = fopen(fileName, 'r', 'ieee-be.l64');
fileHeader = readFileHeader(f);
% if ~strcmp(fileHeader.moduleType, 'WhistlesMoans') || ...
%     ~strcmp(fileHeader.streamName, 'Clicks')
%     disp(sprintf('File %s is not from the Click Detector Detector', fileName));
% end
bad = 0;
iOb = 0;
prevPos = -1;
    iNoise = 0;
while (1)
    iOb = iOb + 1;
    pos = ftell(f);
    if (pos == prevPos)
        disp(sprintf('File stuck at byte %d', pos));
        break;
    end
    prevPos = pos;

    [nextLen, nL] = fread(f, 1, 'int32');
    [nextType, nT] = fread(f, 1, 'int32');
    if (nL == 0 || nT == 0) 
        break;
    end
    fseek(f, -8, 'cof');
    switch nextType
        case -1
            fileHeader = readFileHeader(f);
        case -2
            fileFooter = readFileFooter(f);
        case -3
            moduleHeader = readNoiseModuleHeader(f);
        case -4
            moduleFooter = readModuleFooter(f);
        otherwise
            if (nextType == 1)
                fseek(f, 8, 'cof');
                noise.millis =  fread(f, 1, 'int64');
                noise.date = millisToDateNum(noise.millis);
                noiseDataLen = fread(f, 1, 'int32');
                noise.iChan = fread(f, 1, 'int16');
                noise.nBands = fread(f, 1, 'int16');
                noise.nMeasures = fread(f, 1, 'int16');
                n = fread(f, noise.nBands*noise.nMeasures, 'int16') / 100.;
                noise.noise = reshape(n, noise.nMeasures, noise.nBands);
                if (noise.date > endDate)
                    break;
                end
                if (noise.date >= startDate)
                    iNoise = iNoise + 1;
                    if (iNoise == 1)
                        clear noises;
                        if (preallocate > 0)
                            noises(preallocate) = noise;
                        end
                    end
                    noises(iNoise) = noise;
                end
            else
                fseek(f, nextLen, 'cof');
            end
    end         
end
catch mError
    disp('Error reading file');
end
if preallocate
    if iNoise < preallocate
        noises = noises(1:iNoise);
    end
end
fclose(f);
end
function header = readNoiseModuleHeader(f)
header.length = fread(f, 1, 'int32');
header.identifier = fread(f, 1, 'int32');
header.version = fread(f, 1, 'int32');
header.binaryLength = fread(f, 1, 'int32');
if (header.binaryLength == 0)
    return
end
if (header.binaryLength >= 4)
   header.nBands = fread(f, 1, 'int16'); 
   header.statsTypes = fread(f, 1, 'int16'); 
end
    header.loEdges = fread(f, header.nBands, 'float');
    header.hiEdges = fread(f, header.nBands, 'float');
% fseek(file, header.binaryLength, 'cof');
% fseek(file, 0, 'cof');
end