function ok = writeModuleHeader(fid, version)
if isnumeric(version)
    % assume an empty header with this version number
    headLength = 16;
    ok = true;
    try
        fwrite(fid, headLength, 'int32');
        fwrite(fid, -3, 'int32');
        fwrite(fid, version, 'int32');
        fwrite(fid, 0, 'int32');
    catch
        ok = false;
    end
else
    % assume it's a real header with the correct info in it.
    % header.length = fread(file, 1, 'int32');
    % header.identifier = fread(file, 1, 'int32');
    % header.version = fread(file, 1, 'int32');
    % header.binaryLength = fread(file, 1, 'int32');
    % if readBinary
    %     header.binData = fread(file, header.binaryLength, 'int8');
    ok = true;
    try
        fwrite(fid, version.length, 'int32');
        fwrite(fid, version.identifier, 'int32');
        fwrite(fid, version.version, 'int32');
        fwrite(fid, version.binaryLength, 'int32');
        fwrite(fid, version.binData, 'int8');
    catch
        disp('Error writing module header');
        ok = false;
    end
end
% header.length = fread(file, 1, 'int32');
% header.identifier = fread(file, 1, 'int32');
% header.version = fread(file, 1, 'int32');
% header.binaryLength = fread(file, 1, 'int32');
% fseek(file, header.binaryLength, 'cof');
% fseek(file, 0, 'cof');