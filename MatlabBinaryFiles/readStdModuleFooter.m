function footer = readStdModuleFooter(file)
% reads the module footer information common to all modules.  Differs from
% the legacy code in that it does not read in or skip any information
% specific to a module.  

footer.length = fread(file, 1, 'int32');
footer.identifier = fread(file, 1, 'int32');
footer.binaryLength = fread(file, 1, 'int32');
