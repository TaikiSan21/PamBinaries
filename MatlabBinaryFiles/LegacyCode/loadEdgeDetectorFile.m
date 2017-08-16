function [tones fileHeader fileFooter] = loadEdgeDetectorFile(fileName)
tones = [];

f = fopen(fileName, 'r', 'ieee-be.l64');
fileHeader = readFileHeader(f);
if ~strcmp(fileHeader.moduleType, 'RW Edge Detector') || ...
    ~strcmp(fileHeader.streamName, 'Edges')
    disp(sprintf('File %s is not from the Right whale edge detector', fileName));
end
if strcmp(fileHeader.branch, 'PamBuoy') == 0
end
bad = 0;
iOb = 0;
prevPos = -1;
sampleRate = 2000;
fftLength = 256;
try
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
            moduleheader = readModuleHeader(f);
        case -4
            moduleFooter = readModuleFooter(f);
        case -5
            fseek(f, nextLen, 'cof');
        otherwise
            if (nextType == 0)
                fseek(f, 8, 'cof');
                wsl.millis = fread(f, 1, 'int64');
                wsl.datenum = millisToDateNum(wsl.millis);
                dataLength = fread(f, 1, 'int32');
                wsl.startSample = fread(f, 1, 'int64');
                wsl.channelMap = fread(f, 1, 'int32');
                wsl.type = fread(f, 1, 'int16');
                wsl.signal = fread(f, 1, 'float');
                wsl.noise = fread(f, 1, 'float');
                wsl.nSlices = fread(f, 1, 'int16');
                
                wsl.times = zeros(1, wsl.nSlices);
                wsl.sliceNums = zeros(1, wsl.nSlices);
                wsl.loFreqs = zeros(1, wsl.nSlices);
                wsl.peakFreqs = zeros(1, wsl.nSlices);
                wsl.hiFreqs = zeros(1, wsl.nSlices);
                wsl.peakAmps = zeros(1, wsl.nSlices);
                for i = 1:wsl.nSlices
                    wsl.sliceNums(i) = fread(f, 1, 'int16');
                    wsl.loFreqs(i) = fread(f, 1, 'int16');
                    wsl.peakFreqs(i) = fread(f, 1, 'int16');
                    wsl.hiFreqs(i) = fread(f, 1, 'int16');
                    wsl.peakAmps(i) = fread(f, 1, 'float');
                end
                wsl.loFreqs = wsl.loFreqs*sampleRate/fftLength;
                wsl.peakFreqs = wsl.peakFreqs*sampleRate/fftLength;
                wsl.hiFreqs = wsl.hiFreqs*sampleRate/fftLength;

                tones = [tones wsl];
            end
    end
end
catch mError
    mError
end
fclose(f);