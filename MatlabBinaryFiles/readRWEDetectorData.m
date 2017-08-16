function [data, error] = readRWEDetectorData(fid, fileInfo, data)
% reads binary data stored by the Right Whale Edge Detector.
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

    if (fileInfo.moduleHeader.version==0)
        data.startSample = fread(fid, 1, 'int64');
        data.channelMap = fread(fid, 1, 'int32');
    end
    
    data.type = fread(fid, 1, 'int16');
    data.signal = fread(fid, 1, 'float');
    data.noise = fread(fid, 1, 'float');
    data.nSlices = fread(fid, 1, 'int16');
    
    data.times = zeros(1, data.nSlices);
    data.sliceNums = zeros(1, data.nSlices);
    data.loFreqs = zeros(1, data.nSlices);
    data.peakFreqs = zeros(1, data.nSlices);
    data.hiFreqs = zeros(1, data.nSlices);
    data.peakAmps = zeros(1, data.nSlices);
    for i = 1:data.nSlices
        data.sliceNums(i) = fread(fid, 1, 'int16');
        data.loFreqs(i) = fread(fid, 1, 'int16');
        data.peakFreqs(i) = fread(fid, 1, 'int16');
        data.hiFreqs(i) = fread(fid, 1, 'int16');
        data.peakAmps(i) = fread(fid, 1, 'float');
    end
    
catch mError
    disp(['Error reading ' fileInfo.fileHeader.moduleType '  data object.  Data read:']);
    disp(data);
    disp(getReport(mError));
    error=true;
end
    