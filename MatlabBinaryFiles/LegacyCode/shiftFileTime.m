function shiftFileTime(srcFile, dstFolder, timeShiftMillis)
% Function to shift the time of a PAMguard binary file, but otherwise leave
% it intact. It has to be assumed that all time information is in the
% standard object headers since it won't get altered anywhere else in the
% data
% srcFile - source file
% dstFolder - destination folder. The new file will have a new name since
% it's time will have changed
% timeShiftMillis Time shift in milliseconds. This will be added to all
% times in the data.
if nargin == 0
    srcFile = 'C:\\DecimusTest\\EOM\\binData\\20150122\\Noise_Band_Noise_bands_Noise_bands_20150122_150000.pgdf';
    %     srcFile = 'C:\\DecimusTest\\EOM\\binData\\20150122\\Noise_Band_Noise_bands_Noise_bands_20150122_130000.pgdf';
    dstFolder =  'C:\\DecimusTest\\EOM\\binDataTShift';
    timeShiftMillis = 3600*1000;
end
% will first have to extract the file header.
addDays = timeShiftMillis / (3600*24*1000);
src = fopen(srcFile, 'r', 'ieee-be.l64');
try
    fileHeader = readFileHeader(src, true);
catch
    disp(sprintf('File %s does not contain a valid header and will be skipped', srcFile));
    return;
end
outName = createPamFileName(dstFolder, fileHeader.moduleType, fileHeader.moduleName, ...
    fileHeader.streamName, fileHeader.dataDate+addDays);
% check the output folder exists
a = dir(dstFolder);
if numel(a) == 0
    try
        mkdir(dstFolder);
    catch
        disp(sprintf('Unable to create folder %s', dstFolder));
        return;
    end
end
% now make a file in that folder
dst = fopen(outName, 'w', 'ieee-be.l64');
fileHeader.dataDate = fileHeader.dataDate + addDays;
writeFileHeader(dst, fileHeader);
nobj = 0;

while feof(src) == 0
   % read an object. 
    [nextLen, nL] = fread(src, 1, 'int32');
    [nextType, nT] = fread(src, 1, 'int32');
    if (nL == 0 || nT == 0) 
        break;
    end
    fseek(src, -8, 'cof');
    switch nextType
        case -1
            fileHeader = readFileHeader(src, true);
            fileHeader.dataDate = fileHeader.dataDate + addDays;
            writeFileHeader(dst, fileHeader);
        case -2
            fileFooter = readFileFooter(src);
            fileFooter.dataDate = fileFooter.dataDate + addDays;
            writeFileFooter(dst, fileFooter);
        case -3
            moduleHeader = readModuleHeader(src, true);
            writeModuleHeader(dst, moduleHeader);
        case -4
            moduleFooter = readModuleFooter(src, true);
            writeModuleFooter(dst, moduleFooter);
        otherwise 
            % extract and write data. 
            fseek(src, 8, 'cof');
            fwrite(dst, nextLen, 'int32');
            fwrite(dst, nextType, 'int32');
            millis =  fread(src, 1, 'int64');
            fwrite(dst, millis+timeShiftMillis, 'int64');
            dataLen = fread(src, 1, 'int32');
            fwrite(dst, dataLen, 'int32');
            data = fread(src, dataLen, 'int8');
            fwrite(dst, data, 'int8'); 
            
            nobj = nobj + 1;
    end
    
end

fclose(src);
fclose(dst);
end
