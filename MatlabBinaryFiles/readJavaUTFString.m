function [str len] = readJavaUTFString(file);
% function [str len] = readJavaUTFString(file);
% read a string written in Java UTF-8 format from a file.
% The first two bytes are the length of the string, then 
% it's the string. 
len = fread(file, 1, 'int16');
str = char(fread(file, len, 'uchar')');