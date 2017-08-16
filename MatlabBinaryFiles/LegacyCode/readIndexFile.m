function [fileheader filefooter moduleheader modulefooter] = readIndexFile(fileName)
fileheader  = [];
filefooter = [];
moduleheader = [];
modulefooter = [];
f = fopen(fileName, 'r', 'ieee-be.l64');
while (true)
    [nextLen, nL] = fread(f, 1, 'int32');
    [nextType, nT] = fread(f, 1, 'int32');
    if (nL == 0 || nT == 0) 
        break;
    end
    fseek(f, -8, 'cof');
    switch nextType
        case -1
            fileheader = readFileHeader(f);
        case -2
            filefooter = readFileFooter(f);
        case -3
            moduleheader = readModuleHeader(f);
        case -4
            modulefooter = readModuleFooter(f);
        otherwise
            fseek(f, nextLen, 'cof');
    end         
end
fclose(f);