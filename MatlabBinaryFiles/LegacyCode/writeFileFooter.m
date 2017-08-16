function ok = writeFileFooter(file, footer)
ok = true;
try
    fwrite(file, footer.length, 'int32');
    fwrite(file, footer.identifier, 'int32');
    fwrite(file, footer.nObjects, 'int32');
    fwrite(file, dateNumToMillis(footer.dataDate), 'int64');
    fwrite(file, dateNumToMillis(footer.analysisDate), 'int64');
    fwrite(file, footer.endSample, 'int64');
    fwrite(file, footer.fileLength, 'int64');
    fwrite(file, footer.endReason, 'int32');
catch
    ok = false;
end