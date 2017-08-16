function [dbhts fileHeader fileFooter] = loadDBHTFile(fileName, arg2, arg3, arg4)
dbhts = [];
fileHeader = [];
fileFooter = [];
if nargin == 0
    error('A file name is required by loadDBHTfile')
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
    iDBHT = 0;
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
                if (nextType == 1)
                    fseek(f, 8, 'cof');
                    dbht.millis =  fread(f, 1, 'int64');
                    dbht.date = millistoDateNum(dbht.millis);
                    dbht.startSample = fread(f, 1, 'int64');
                    fread(f, 1, 'int32');
                    dbht.channelMap = fread(f, 1, 'int32');
                    dbht.rms = fread(f, 1, 'int16')/100;
                    dbht.zeroPeak = fread(f, 1, 'int16')/100;
                    dbht.peakPeak = fread(f, 1, 'int16')/100;
                    if (dbht.date > endDate)
                        break;
                    end
                    if (dbht.date >= startDate && dbht.date <= endDate)
                        iDBHT = iDBHT + 1;
                        if (iDBHT == 1)
                            clear dbhts;
                            if (preallocate > 0)
                                dbhts(preallocate) = dbht;
                            end
                        end
                        dbhts(iDBHT) = dbht;
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
    if iDBHT < preallocate
        dbhts = dbhts(1:iClick);
    end
end
try
    fclose(f);
catch
    disp('Bad file close')
end