function [data, error] = readNoiseBandData(fid, fileInfo, data)
% reads binary data stored by the Noise Band Monitor.
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
    
    % read module specific data
    dataLength = fread(fid, 1, 'int32');
    if (dataLength==0)
        return;
    end

    if (fileInfo.moduleHeader.version<=2)
        data.startSample = fread(fid, 1, 'int64');
        data.channelMap = fread(fid, 1, 'int32');
    end
    
    data.rms = fread(fid, 1, 'int16')/100.;
    data.zeroPeak = fread(fid, 1, 'int16')/100.;
    data.peakPeak = fread(fid, 1, 'int16')/100.;
    
    if (fileInfo.moduleHeader.version>=2)    
        data.sel = fread(fid, 1, 'int16')/100.;
        data.selSecs = fread(fid, 1, 'int16');
    end
    
catch mError
    disp(['Error reading ' fileInfo.fileHeader.moduleType '  data object.  Data read:']);
    disp(data);
    disp(getReport(mError));
    error=true;
end
    