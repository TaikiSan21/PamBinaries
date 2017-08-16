function data = readPamData(fid, fileInfo)
% Reads in the object data that is common to all modules.  This reads up to
% (but not including) the object binary length, and then calls a function
% to read the module-specific data.

% Inputs:
%   fid = file identifier
%   fileInfo = structure holding the file header, module header, a handle
%   to the function to read module-specific data, etc.
%
% Output:
%   data = structure containing data from a single object
%


% set constants to match flag bitmap constants in class
% DataUnitBaseData.java.  The following contstants match header version 4.
TIMEMILLIS          = hex2dec('1');
TIMENANOS           = hex2dec('2');
CHANNELMAP          = hex2dec('4');
UID                 = hex2dec('8');
STARTSAMPLE         = hex2dec('10');
SAMPLEDURATION      = hex2dec('20');
FREQUENCYLIMITS     = hex2dec('40');
MILLISDURATION      = hex2dec('80');
TIMEDELAYSSECS      = hex2dec('100');

% initialize a new variable to hold the data
data=[];
data.flagBitmap=0;

% calculate where the next object starts, in case there is an error trying
% to read this one
objectLen = fread(fid, 1, 'int32');
curObj = ftell(fid);
nextObj = curObj + objectLen;

% first thing to check is that this is really the type of object we think
% it should be, based on the file header.  If not, warn the user, move the
% pointer to the next object, and exit
data.identifier = fread(fid, 1, 'int32');
if (any(data.identifier==fileInfo.objectType))
    % do nothing here - couldn't figure out a clean way of checking if
    % number wasn't in array
else
    disp(['Error - Object Identifier does not match ' fileInfo.fileHeader.moduleType ' type.  Aborting data read.']);
    fseek(fid, nextObj, 'bof');
    return;
end
    
% read the data, starting with the standard data that every data unit has
version=fileInfo.fileHeader.fileFormat;
try
    data.millis = fread(fid, 1, 'int64');
    
    if (version >=3)
        data.flagBitmap = fread(fid,1,'int16');
    end
    
    if (version == 2 || (bitand(data.flagBitmap, TIMENANOS)~=0) )
        data.timeNanos = fread(fid,1,'int64');
    end
    
    if (version == 2 || (bitand(data.flagBitmap, CHANNELMAP)~=0) )
        data.channelMap = fread(fid,1,'int32');
    end
    
    if (bitand(data.flagBitmap, UID)==UID)
        data.UID = fread(fid,1,'int64');
    end
    
    if (bitand(data.flagBitmap, STARTSAMPLE)~=0)
        data.startSample = fread(fid,1,'int64');
    end
    
    if (bitand(data.flagBitmap, SAMPLEDURATION)~=0)
        data.sampleDuration = fread(fid,1,'int32');
    end
    
    if (bitand(data.flagBitmap, FREQUENCYLIMITS)~=0)
        minFreq = fread(fid,1,'float');
        maxFreq = fread(fid,1,'float');
        data.freqLimits = [minFreq maxFreq];
    end
    
    if (bitand(data.flagBitmap, MILLISDURATION)~=0)
        data.millisDuration = fread(fid,1,'float');
    end
    
    if (bitand(data.flagBitmap, TIMEDELAYSSECS)~=0)
        data.numTimeDelays = fread(fid,1,'int16');
        td=zeros(data.numTimeDelays);
        for i = 1:data.numTimeDelays
            td(i)=fread(fid,1,'float');
        end
        data.timeDelays=td;
    end
    
    % set date, to maintain backwards compatibility
    data.date = millisToDateNum(data.millis);
    
    % now read the module-specific data
    if(isa(fileInfo.readModuleData,'function_handle'))
        [data, error] = fileInfo.readModuleData(fid, fileInfo, data);
        if (error)
            disp(['Error - cannot retrieve ' fileInfo.fileHeader.moduleType ' data properly.']);
            fseek(fid, nextObj, 'bof');
            return;
        end
    end
    
    
catch mError
    disp('Error loading object data');
    disp(data);
    disp(getReport(mError));
    fseek(fid, nextObj, 'bof');
end