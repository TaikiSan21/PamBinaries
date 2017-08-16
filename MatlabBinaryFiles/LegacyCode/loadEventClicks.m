function evClicks = loadEventClicks(dataDir, fileNames, clickNos, eventDateRange)
% function to load all clicks from an event even if spread over multiple files. 
% inputs are root data folder, cell array of binary file names and click
% numbers
evClicks = [];
if nargin < 4
    eventDateRange = [0 datenum(2100,1,1)];
end

% first a bit of mucking about to get a directory listing of all 
% files. These are held globally so that they can be reused in 
% subsequent calls. 
global openClickDir openFileList openFileNames;
if ~isempty(openClickDir)
    if strcmp(openClickDir, dataDir) == 0
        openClickDir = [];
    end
end
fileNames = deblank(fileNames);
if isempty(openClickDir)
    % remake the file list. 
    if numel(fileNames) == 0
        error('null file list in loadEventClicks')
    end
    file1 = fileNames{1};
    filePattern = [file1(1:(end-21)) '*.pgdf'];
    openFileList = dirsub(dataDir, filePattern);
    openClickDir = dataDir;
    % now make a list of the file names only
    for i = numel(openFileList):-1:1
        [path name ext] = fileparts(openFileList(i).name);
        openFileNames{i} = [name ext];
    end
end
% should now be able to find files pretty efficiently !
% unique file list
unFiles = unique(fileNames);
for i = 1:numel(unFiles) % loop over the different files
    % list of clicks in a particular file
    fileClickNos = clickNos(find(strcmp(fileNames, unFiles{i})));
    % now find the file !
    fileIndex = find(strcmp(openFileNames, unFiles{i}));
    if numel(fileIndex) ~= 1
        error(sprintf('Unable to find click file %s', unFiles{i}));
    end
    % put start date 0 to always read from start of file - makes indexing
    % easier. 
    fileClicks = loadClickFile(openFileList(fileIndex).name, ...
         0, eventDateRange(2), 5000);
    % click numbers are zero indexed from the database but 1 indexed from file.
    evClicks = [evClicks fileClicks(fileClickNos+1)];
end

