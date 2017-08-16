function [fileName indexName] = createPamFileName(path, moduleType, moduleName, streamName, date)
% function fileName = createPamFileName(path, moduleType, moduleName, streamName, date) 
% create a PAMGUARD binary file name compatible with the standard PAMGUARD
% binary storage system based on the path, moduleType, moduleName,
% streamName and date
if isempty(path) || length(path) == 0
    path = '.\';
elseif path(end) ~= '\'
    path = [path '\'];
end
nameBase = sprintf('%s%s_%s_%s_%s', path, moduleType, moduleName, ...
    streamName,  datestr(date, 'yyyymmdd_HHMMSS'));
nameBase = strrep(nameBase, ' ', '_');
fileName = [nameBase '.pgdf'];
if nargout >= 2
    indexName = [nameBase '.pgdx'];
end