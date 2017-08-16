function ok = writeFileHeader(file, header)
ok = 1;
try
    fwrite(file, header.length, 'int32');
    fwrite(file, header.identifier, 'int32');
    fwrite(file, header.fileFormat, 'int32');
    fwrite(file, real('PAMGUARDDATA'), 'uchar');
    writeJavaUTFString(file, header.version);
    writeJavaUTFString(file, header.branch);
    fwrite(file, dateNumToMillis(header.dataDate), 'int64');
    fwrite(file, dateNumToMillis(header.analysisDate), 'int64');
    fwrite(file, header.startSample, 'int64');
    writeJavaUTFString(file, header.moduleType);
    writeJavaUTFString(file, header.moduleName);
    writeJavaUTFString(file, header.streamName);
    fwrite(file, header.extraInfoLen, 'int32');
    if (header.extraInfoLen > 0)
    fwrite(file, header.extraInfo, 'int8');
    end
    % fseek(file, header.extraInfoLen, 'cof');
catch
    ok = false
end
