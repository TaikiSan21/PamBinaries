function [daqStat fileHeader fileFooter moduleheader modulefooter] = loadDaqStatusFile(fileName, arg2, arg3, arg4)
daqStat = [];
fileHeader = [];
fileFooter = [];
moduleheader = [];
modulefooter = [];
if nargin == 0
    error('A file name is required by loadDaqStatusFile')
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
    % conversion constants to convert values
    a = 127.*2./log(32767);
    b = -127;
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
                moduleheader = readDaqStatModuleHeader(f);
                if (moduleheader.fftLength == 0) 
                    return
                end
            case -4
                moduleFooter = readModuleFooter(f);
            otherwise
                if (nextType == 1)
                    fseek(f, 8, 'cof');
                    ltsa.millis =  fread(f, 1, 'int64');
                    ltsa.date = millisToDateNum(ltsa.millis);
                    ltsa.startSample = fread(f, 1, 'int64');
                    fread(f, 1, 'int32');
                    ltsa.channelMap = fread(f, 1, 'int32');
                    ltsa.endMillis = fread(f, 1, 'int64');
                    ltsa.endDate = millisToDateNum(ltsa.endMillis);
                    ltsa.nFFT = fread(f, 1, 'int32');
                    ltsa.maxVal = fread(f, 1, 'single');
                    ltsa.byteData = fread(f, moduleheader.fftLength/2, 'int8');
                    % convert format...
                    ltsa.data = exp((ltsa.byteData-b)/a)*ltsa.maxVal / 32767;
                    if (ltsa.date > endDate)
                        break;
                    end
                    if (ltsa.date >= startDate && ltsa.date <= endDate)
                        iDBHT = iDBHT + 1;
                        if (iDBHT == 1)
                            clear ltsas;
                            if (preallocate > 0)
                                ltsas(preallocate) = ltsa;
                            end
                        end
                        ltsas(iDBHT) = ltsa;
                    end
                else
                    fseek(f, nextLen, 'cof');
                end
        end
    end
catch mError
    disp('Error reading file');
    disp(mError.message);
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
end

function header = readDaqStatModuleHeader(file);
header.length = fread(file, 1, 'int32');
header.identifier = fread(file, 1, 'int32');
header.version = fread(file, 1, 'int32');
header.binaryLength = fread(file, 1, 'int32');
header.deviceName = readJavaUTFString(file);
header.sampleRate = fread(file, 1, 'float32');
header.fftLength = fread(file, 1, 'int32');
header.fftHop = fread(file, 1, 'int32');
header.intervalSeconds = fread(file, 1, 'int32');
% fseek(file, header.binaryLength, 'cof');
end
    