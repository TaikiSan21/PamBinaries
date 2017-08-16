function [clicks fileHeader fileFooter] = loadClickFile(fileName, arg2, arg3, arg4)
clicks = [];
fileHeader = [];
fileFooter = [];
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
iClick = 0;
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
    fileClickNo = 0;
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
            otherwise
                if (nextType == 1000)
                    fseek(f, 8, 'cof');
                    fileClickNo = fileClickNo + 1; % go for one indexing.
                    click.clickNo = fileClickNo;
                    click.millis =  fread(f, 1, 'int64');
                    click.date = millisToDateNum(click.millis);
                    clickDataLen = fread(f, 1, 'int32');
                    click.startSample = fread(f, 1, 'int64');
                    click.channelMap = fread(f, 1, 'int32');
                    click.triggerMap = fread(f, 1, 'int32');
                    click.type = fread(f, 1, 'int16');
                    if (moduleheader.version >= 3)
                        click.flags = fread(f, 1, 'int32');
                    else
                        click.flags = 0;
                    end
                    nDelays = fread(f, 1, 'int16');
                    if (nDelays)
                        click.delays = fread(f, nDelays, 'float');
                    end
                    nAngles = fread(f, 1, 'int16');
                    if (nAngles)
                        click.angles = fread(f, nAngles, 'float');
                    end
                    if (moduleheader.version >= 3)
                        nAngleErrors = fread(f, 1, 'int16');
                        click.angleErrors = fread(f, nAngleErrors, 'float');
                    else
                        click.angleErrors = [];
                    end
                    click.duration = fread(f, 1, 'int16');
                    click.nChan = countChannels(click.channelMap);
                    maxVal = fread(f, 1, 'float');
                    click.wave = fread(f, [click.duration,click.nChan], 'int8') * maxVal / 127;
                    if (click.date > endDate)
                        break;
                    end
                    if (click.date >= startDate)
                        iClick = iClick + 1;
                        if (iClick == 1)
                            clear clicks;
                            if (preallocate > 0)
                                clicks(preallocate) = click;
                            end
                        end
                        clicks(iClick) = click;
%                         if mod(iClick,1000) == 0
%                             fprintf('%d clicks loaded from file\n', iClick);
%                         end
                    end
                else
                    fseek(f, nextLen, 'cof');
                end
        end
    end
catch mError
    disp(['Error reading file: ' mError.message]);
end
if preallocate & numel(clicks) > 0
    if iClick < preallocate
        clicks = clicks(1:iClick);
    end
end
fclose(f);