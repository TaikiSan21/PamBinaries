function [data, error] = readNoiseMonData(fid, fileInfo, data)
% reads binary data stored by the Noise Monitor.
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

    data.iChan = fread(fid, 1, 'int16');
    data.nBands = fread(fid, 1, 'int16');
    
    if (fileInfo.moduleHeader.version>=1)    
        data.nMeasures = fread(fid, 1, 'int16');
    else
        data.nMeasures = 4;
    end
    
    if (fileInfo.moduleHeader.version<=1)
        n = fread(fid, data.nBands*data.nMeasures, 'float');
    else
        n = fread(fid, data.nBands*data.nMeasures, 'int16') / 100.;
    end
    data.noise = reshape(n, data.nMeasures, data.nBands);
    
    
    
catch mError
    disp(['Error reading ' fileInfo.fileHeader.moduleType '  data object.  Data read:']);
    disp(data);
    disp(getReport(mError));
    error=true;
end
    