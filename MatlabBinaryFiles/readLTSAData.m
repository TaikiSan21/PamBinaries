function [data, error] = readLTSAData(fid, fileInfo, data)
% reads binary data stored by the LTSA module.
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
a = 127.*2./log(32767);
b = -127;

try
    
    % read LTSA specific data
    dataLength = fread(fid, 1, 'int32');
    if (dataLength==0)
        return;
    end

    if (fileInfo.moduleHeader.version<=1)
        data.startSample = fread(fid, 1, 'int64');
    end

    if (fileInfo.moduleHeader.version==0)
        data.duration = fread(fid, 1, 'int64');
    end
    
    if (fileInfo.moduleHeader.version<=1)
        data.channelMap = fread(fid, 1, 'int32');
    end
    
    data.endMillis = fread(fid, 1, 'int64');
    data.endDate = millisToDateNum(data.endMillis);
    data.nFFT = fread(fid, 1, 'int32');
    data.maxVal = fread(fid, 1, 'float');
    
    % version 0 scaled the data linearly to 16 bit
    if (fileInfo.moduleHeader.version==0)    
        data.byteData = fread(fid, fileInfo.moduleHeader.fftLength/2, 'int16');
        data.data = data.byteData / 32767. * data.maxVal;
        
    % after version 0, the data was first scaled to 16 bit and then
    % converted to a log so that it could be saved as an 8 bit.
    else
        data.byteData = fread(fid, fileInfo.moduleHeader.fftLength/2, 'int8');
        data.data = exp((data.byteData-b)/a)*data.maxVal / 32767;
    end
    
catch mError
    disp(['Error reading ' fileInfo.fileHeader.moduleType '  data object.  Data read:']);
    disp(data);
    disp(getReport(mError));
    error=true;
end

    