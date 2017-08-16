function footer = readClickFooter(file)
% reads module footer information for the Click Detector module.  Note that
% sometimes there is no additional footer information, so check first
% whether or not the binaryLength variable is 0.

footer=readStdModuleFooter(file);
if (footer.binaryLength ~= 0)
    footer.typesCountLength = fread(file, 1, 'int16');
    footer.typesCount = fread(file, footer.typesCountLength, 'int32');
end
