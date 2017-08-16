function shiftFileTimes(srcFolder, dstFolder, timeShiftMillis)
% shift a load of binary file data times by the given amount of
% milliseconds. Can be used when time zone disasters occurr !!
% srcFolder - source folder
% dstFolder - destination folder
% timeShiftNillis - time shift in milliseconds. 

% this function should now be available on SVN
% will automatically scan subfolders of dstFolder and add the same sub
% folder structure to dstFolder.

if nargin == 0
    srcFolder = 'C:\DecimusTest\EOM\binData';
    dstFolder = 'C:\DecimusTest\EOM\binDataTShift';
    timeShiftMillis = 3600000;
end
mkdir(dstFolder);
allSrcFiles = dirsub(srcFolder, '*.pgdf');
for i = 1:numel(allSrcFiles)
    aFile = allSrcFiles(i).name;
    aPath = fileparts(aFile);
    pathEnd = aPath(length(srcFolder)+1:end);
    aDest = [dstFolder pathEnd];
    shiftFileTime(aFile, aDest, timeShiftMillis);
end