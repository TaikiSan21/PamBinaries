function [moduleNameData, settingsList, header, footer] = readSettingsFile(fileName)
f = fopen(fileName, 'r', 'ieee-be.l64');
moduleNameData = {};
settingsList = {};
header = readFileHeader(f);
if strcmp(header.moduleType, 'SettingsStore') == 0
    error('opened file is not a pamguard settings store')
end
while (1)
    nextObjectStart = ftell(f);
    nextObjectSize = fread(f, 1, 'int32');
    nextObjectId = fread(f, 1, 'int32');
    if (nextObjectId < 0)
        break;
    end
    nextBinaryLen = fread(f, 1, 'int32');
    if nextObjectId == 1
        moduleNameData{end+1} = readModuleNameData(f, nextBinaryLen);
    elseif nextObjectId == 2
        settingsList{end+1} = readSettingsData(f, nextBinaryLen);
    else
        disp(sprintf('unknown object in file with id %d', nextObjectId));
        fseek(f, nextBinaryLen, 'cof');
    end
end
fseek(f, nextObjectStart, 'bof');
footer = readFooter(f);
fclose(f);
end

function moduleNameData = readModuleNameData(f, size)
moduleNameData.class = readJavaUTFString(f);
moduleNameData.moduleType = readJavaUTFString(f);
moduleNameData.moduleName = readJavaUTFString(f);
%fseek(f, size, 'cof');
end

function unitData = readSettingsData(f, size)
unitData.unitType = readJavaUTFString(f);
unitData.unitName = readJavaUTFString(f);
unitData.versionNo = fread(f, 1, 'int64');
unitData.dataSize = fread(f, 1, 'int32');
unitData.data = fread(f, unitData.dataSize, 'uchar')
%fseek(f, unitData.dataSize, 'cof');
end