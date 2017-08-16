function [data, error] = readClickData(fid, fileInfo, data)
% reads binary data stored by the Click Detector.
%
% Inputs:
%   fid = file identifier
%   fileInfo = structure holding the file header, module header, a handle
%   data = a structure containing the standard data
%
% Output:
%   data = structure containing data from a single object
%

% initialize variables
error=false;

try
    
    % read click detector specific data
    dataLength = fread(fid, 1, 'int32');
    if (dataLength==0)
        return;
    end
    
    if (fileInfo.moduleHeader.version<=3)
        data.startSample = fread(fid, 1, 'int64');
        data.channelMap = fread(fid, 1, 'int32');
    end
    
    data.triggerMap = fread(fid, 1, 'int32');
    data.type = fread(fid, 1, 'int16');
    if (fileInfo.moduleHeader.version >= 2)
        data.flags = fread(fid, 1, 'int32');
    else
        data.flags = 0;
    end
    
    if (fileInfo.moduleHeader.version <= 3)
        nDelays = fread(fid, 1, 'int16');
        if (nDelays)
            data.delays = fread(fid, nDelays, 'float');
        end
    end
    
    nAngles = fread(fid, 1, 'int16');
    if (nAngles)
        data.angles = fread(fid, nAngles, 'float');
    end
    
    if (fileInfo.moduleHeader.version >= 3)
        nAngleErrors = fread(fid, 1, 'int16');
        data.angleErrors = fread(fid, nAngleErrors, 'float');
    else
        data.angleErrors = [];
    end
    
    if (fileInfo.moduleHeader.version <= 3)    
        data.duration = fread(fid, 1, 'int16');
    else
        data.duration = data.sampleDuration;    % duplicate the value to maintain backwards compatibility
    end
    
    data.nChan = countChannels(data.channelMap);
    maxVal = fread(fid, 1, 'float');
    data.wave = fread(fid, [data.duration,data.nChan], 'int8') * maxVal / 127;
    
catch mError
    disp(['Error reading ' fileInfo.fileHeader.moduleType '  data object.  Data read:']);
    disp(data);
    disp(getReport(mError));
    error=true;
end
