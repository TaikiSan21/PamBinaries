function map = createDataMap(moduleType, moduleName, streamName, binaryStore, verbose)
% function map = createDataMap(ModuleType, ModuleName, streamName,
% binaryStore, verbose)
% create a data map analagous to the one in PAMGUARD which will speed up
% finding data from specific time periods. Wherever possible use the index
% files but when they don't exist, use the binary files. From each extract
% the header and footer and store in an array of data objects. 
if nargin < 5
    verbose = 0;
end
map = [];
mask = sprintf('%s_%s_%s*.pgdf', moduleType, moduleName, ...
    streamName);
mask = strrep(mask, ' ', '_');
files = dirsub(binaryStore, mask);
for i = 1:length(files)
   dFile =  files(i).name;
   xFile = strrep(dFile, 'pgdf', 'pgdx');
   xD = dir(xFile);
   if length(xD)
       indFile = xFile;
   else
       indFile = dFile;
   end
   [fileheader filefooter moduleheader modulefooter] = readIndexFile(indFile);
   mP.dataFile = dFile;
   mP.indexFile = xFile;
   mP.fileHeader = fileheader;
   mP.fileFooter = filefooter;
   mP.moduleHeader = moduleheader;
   mP.moduleFooter = modulefooter;
   mP.startDate = fileheader.dataDate;
   mP.endDate = 0;
   if ~isempty(filefooter)
       mP.endDate = filefooter.dataDate;
   end
   if i == 1
       clear map
       map(length(files)) = mP;
   end
   map(i) = mP;
   if (verbose & mod(i,1000) == 0)
       disp(sprintf('Reading map point %d of %d = %3.1f%%', i, length(files), i/length(files)*100));
   end
end