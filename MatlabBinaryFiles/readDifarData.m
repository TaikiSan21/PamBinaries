function [data, error] = readDifarData(fid, fileInfo, data)
% reads binary data stored by the Difar module.
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
    
    % read Difar specific data
    dataLength = fread(fid, 1, 'int32');
    if (dataLength==0)
        return;
    end
    
    if (fileInfo.moduleHeader.version <= 1)
        data.startSample = fread(fid, 1, 'int64');
    end

    data.clipStart = fread(fid, 1, 'int64');
    
    if (fileInfo.moduleHeader.version <= 1)
        data.channelMap = fread(fid, 1, 'int32');
    end
    
    data.displaySampleRate = fread(fid, 1, 'float');
    data.demuxedLength = fread(fid, 1, 'int32');
    
    if (fileInfo.moduleHeader.version <= 1)
        minFreq = fread(fid,1,'float');
        maxFreq = fread(fid,1,'float');
        data.freqLimits = [minFreq maxFreq];
    end
    
    data.amplitude = fread(fid, 1, 'float');
    data.gain = fread(fid, 1, 'float');
    data.selAngle = fread(fid, 1, 'float');
    data.selFreq = fread(fid, 1, 'float');
    data.speciesCode = readJavaUTFString(fid);

    if (fileInfo.moduleHeader.version >= 1)
        data.trackedGroup = readJavaUTFString(fid);
    end
    
    data.maxVal = fread(fid, 1, 'float');
    data.demuxData = fread(fid, [data.demuxedLength,3], 'int16') * data.maxVal / 32767;
    
    data.numMatches = fread(fid, 1, 'int16');
    if (data.numMatches > 0)
        data.latitude = fread(fid, 1, 'float');
        data.longitude = fread(fid, 1, 'float');
        
        if (fileInfo.moduleHeader.version >= 1)
            errorX = fread(fid,1,'float');
            errorY = fread(fid,1,'float');
            data.errors = [errorX errorY 0];
        end
        
        for i=1:data.numMatches-1
            data.matchChan(i) = fread(fid, 1, 'int16');
            data.matchTime(i) = fread(fid, 1, 'int64');
        end
    end
        
    
catch mError
    disp(['Error reading ' fileInfo.fileHeader.moduleType '  data object.  Data read:']);
    disp(data);
    disp(getReport(mError));
    error=true;
end
