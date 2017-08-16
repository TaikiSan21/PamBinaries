function header = createPamHeader(moduleType, moduleName, streamName, dataDate)

header.identifier = -1;
header.fileFormat = 1;
header.pamguard = 'PAMGUARDDATA';
header.version = 1.0;
header.branch = 'D3';
header.dataDate = dataDate;
header.analysisDate = now();
header.startSample = 0;
header.moduleType = moduleType;
header.moduleName = moduleName;
header.streamName = streamName;
header.extraInfoLen = 0;
header.length =  12 + length(header.pamguard) + ...
       length(header.version) + ...
		length(header.branch) + 8 + 8 + 8 + ...% three 8's for the three times
		2 + length(header.moduleType) + ...
		2 + length(header.moduleName) + ... 
		2 + length(header.streamName) + ...
		4;