function [data, error] = readClipData(fid, fileInfo, data)
% reads binary data stored by the Clip Generator module.
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
    
    % read Clip Generator specific data
    dataLength = fread(fid, 1, 'int32');
    if (dataLength==0)
        return;
    end
    
    if (fileInfo.moduleHeader.version <= 1)
        data.startSample = fread(fid, 1, 'int64');
        data.channelMap = fread(fid, 1, 'int32');
    end
    
    data.triggerMillis = fread(fid, 1, 'int64');
    
    if (fileInfo.moduleHeader.version <= 1)
        data.sampleDuration = fread(fid, 1, 'int32');
    end
    
    data.filename = readJavaUTFString(fid);
    data.triggerName = readJavaUTFString(fid);
    
    % check if the object type = 2.  If it is, there must be wav data at
    % the end of this object
    if (data.identifier==2)
        data.nChan = fread(fid, 1, 'int16');
        data.nSamps = fread(fid, 1, 'int32');
        data.scale = 1/fread(fid, 1, 'float');
       data.wave = fread(fid, [data.nSamps,data.nChan], 'int8') * data.scale;
    end
    
catch mError
    disp(['Error reading ' fileInfo.fileHeader.moduleType '  data object.  Data read:']);
    disp(data);
    disp(getReport(mError));
    error=true;
end
