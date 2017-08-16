function [data, error] = readAISData(fid, fileInfo, data)
% reads binary data stored by the AIS Processing module.
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
    
    % read AIS Processing Module specific data
    dataLength = fread(fid, 1, 'int32');
    if (dataLength==0)
        return;
    end

    data.mmsiNumber = fread(fid, 1, 'int32');
    data.fillBits = fread(fid, 1, 'int16');
    [strVal, strLen] = readJavaUTFString(fid);
    data.charData = strVal;
    [strVal, strLen] = readJavaUTFString(fid);
    data.aisChannel = strVal;
    
    
catch mError
    disp('Error reading AIS data object.  Data read:');
    disp(data);
    disp(getReport(mError));
    error=true;
end
