function footer = readModuleFooter(file, readBinary)
if nargin < 2
    readBinary = false;
end
footer.length = fread(file, 1, 'int32');
footer.identifier = fread(file, 1, 'int32');
footer.binaryLength = fread(file, 1, 'int32');
if readBinary
    footer.binData = fread(file, footer.binaryLength, 'int8');
else
    footer.binData = [];
    fseek(file, footer.binaryLength, 'cof');
end