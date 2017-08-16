function [tones fileHeader fileFooter moduleHeader moduleFooter] = loadWhistleFile(fileName)
tones = [];
moduleHeader = [];
moduleFooter = [];
fileFooter = [];
fileHeader = [];
if nargin == 0
    fileName = 'C:\WhistleClassifier\48kBS\19900101\WhistlesMoans_Whistle_and_Moan_Detector_Contours_19900101_001000.pgdf'
end

f = fopen(fileName, 'r', 'ieee-be.l64');
fileHeader = readFileHeader(f);
if ~strcmp(fileHeader.moduleType, 'WhistlesMoans') || ...
        ~strcmp(fileHeader.streamName, 'Contours')
    disp(sprintf('File %s is not from the whistle and moan Detector', fileName));
end
if strcmp(fileHeader.branch, 'PamBuoy') == 0
end
bad = 0;
iOb = 0;
prevPos = -1;
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
                moduleHeader = readWhistleModuleHeader(f);
            case -4
                moduleFooter = readModuleFooter(f);
            otherwise
                if (nextType == 2000)
                    fseek(f, 8, 'cof');
                    wsl.millis = fread(f, 1, 'int64');
                    wsl.date = millisToDateNum(wsl.millis);
                    dataLength = fread(f, 1, 'int32');
                    wsl.startSample = fread(f, 1, 'int64');
                    wsl.channelMap = fread(f, 1, 'int32');
                    wsl.nSlices = fread(f, 1, 'int16');
                    if (moduleHeader.version >= 1)
                        wsl.amplitude = fread(f, 1, 'int16') / 100;
                        wsl.nDelays = fread(f, 1, 'int8');
                        wsl.delays = fread(f, wsl.nDelays, 'int16'); % need to scale this still !!!!
                        if ~isempty(moduleHeader)
                            wsl.delays = wsl.delays / moduleHeader.delayScale;
                        end
                    end
                    wsl.sliceData = [];
                    wsl.contour = zeros(1,wsl.nSlices);
                    wsl.contWidth = zeros(1,wsl.nSlices);
                    for i = 1:wsl.nSlices
                        aSlice.sliceNumber = fread(f, 1, 'int32');
                        aSlice.nPeaks = fread(f, 1, 'int8');
                        aSlice.peakData = zeros(4, aSlice.nPeaks);
                        for p = 1:aSlice.nPeaks
                            sss = fread(f, 4, 'int16');
                            aSlice.peakData(:,p) = sss;
                        end
                        wsl.sliceData{i} = aSlice;
                        wsl.contour(i) = aSlice.peakData(2,1);
                        wsl.contWidth(i) = aSlice.peakData(3,1) - aSlice.peakData(1,1) + 1;
                    end
                    wsl.meanWidth = mean(wsl.contWidth);

                    tones = [tones wsl];
                else
                    fseek(f, nextLen-4, 'cof');
                end
        end
    end
catch mError
end
fclose(f);