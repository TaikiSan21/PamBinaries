function [data, error] = readWMDData(fid, fileInfo, data)
% reads binary data stored by the Whistle & Moan Detector.
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
    
   % read WMD specific data
    dataLength = fread(fid, 1, 'int32');
    if (dataLength==0)
        return;
    end
    
    if (fileInfo.moduleHeader.version<=1)
        data.startSample = fread(fid, 1, 'int64');
        data.channelMap = fread(fid, 1, 'int32');
    end
    
    data.nSlices = fread(fid, 1, 'int16');
    
    if (fileInfo.moduleHeader.version >= 1)
        data.amplitude = fread(fid, 1, 'int16') / 100;
    end
    
    if (fileInfo.moduleHeader.version == 1)
        data.nDelays = fread(fid, 1, 'int8');
        data.delays = fread(fid, data.nDelays, 'int16'); % need to scale this still !!!!
        if ~isempty(moduleHeader)
            data.delays = data.delays / moduleHeader.delayScale;
        end
    end
    
    data.sliceData = [];
    data.contour = zeros(1,data.nSlices);
    data.contWidth = zeros(1,data.nSlices);
    for i = 1:data.nSlices
        aSlice.sliceNumber = fread(fid, 1, 'int32');
        aSlice.nPeaks = fread(fid, 1, 'int8');
        aSlice.peakData = zeros(4, aSlice.nPeaks);
        for p = 1:aSlice.nPeaks
            sss = fread(fid, 4, 'int16');
            aSlice.peakData(:,p) = sss;
        end
        data.sliceData{i} = aSlice;
        data.contour(i) = aSlice.peakData(2,1);
        data.contWidth(i) = aSlice.peakData(3,1) - aSlice.peakData(1,1) + 1;
    end
    data.meanWidth = mean(data.contWidth);
    
catch mError
    disp(['Error reading ' fileInfo.fileHeader.moduleType '  data object.  Data read:']);
    disp(data);
    disp(getReport(mError));
    error=true;
end
