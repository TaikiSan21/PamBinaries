function ok = writeJavaUTFString(file, str, len);
% function ok = writeJavaUTFString(file, str, len);
% write a string in UTF-8 format to a file.
% The first two bytes are the length of the string, then
% it's the string.
if nargin < 3
    len = length(str);
end
ok = true;
try
    fwrite(file, len, 'int16');
    fwrite(file, real(str), 'uchar');
catch
    ok = 0;
end