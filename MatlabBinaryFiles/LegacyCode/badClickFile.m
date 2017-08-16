fileName = 'F:\OP\20100804\Click_Detector_Click_Detector_Clicks_20100804_021044.pgdf'
% fileName = 'F:\OP\20100809\Click_Detector_Click_Detector_Clicks_20100809_210000.pgdf'
% fileName = 'F:\OP\20100823\Click_Detector_Click_Detector_Clicks_20100823_044233.pgdf'
fileName = 'C:\PamguardTest\binaryData\Laptop48kHz\20100401\Click_Detector_Click_Detector_Clicks_20100401_161300.pgdf'

fileDir = 'C:\PamguardTest\binaryData\Laptop48kHz\20100824\'
fileDir = 'F:\OP\20100824\'
dr = dir([fileDir '*.pgdf'])
for iF = 1:length(dr)
    fileName = [fileDir dr(iF).name]


    fileSize = dr(iF).bytes;
    if (fileSize < 8) 
        continue;
    end
    file = fopen(fileName, 'r', 'ieee-be.l64');
    fileHeader = readFileHeader(file);
    moduleHeader = readModuleHeader(file)
    bad = 0;
    iOb = 0;
    while (1)
        iOb = iOb + 1;
        pos = ftell(file);
        if (pos >= fileSize)
            break;
        end
        if (feof(file))
            break;
        end
        objectLength = fread(file, 1, 'int32');;
        objectId = fread(file, 1, 'int32');
        if (objectId == 1000)
            fseek(file, objectLength-8, 'cof');
        else
            disp(sprintf('Obj %d, Pos %d/%d, Object type = %d, length = %d', ...
                iOb, pos, fileSize, objectId, objectLength));
            fseek(file, -8, 'cof');
            data = fread(file, objectLength, 'int8')'
        end
        if (objectId == -2)
            break;
        end
        if (objectLength <= 0)
            disp('Zero length object')
            break;
        end

    end

    fclose(file);
end