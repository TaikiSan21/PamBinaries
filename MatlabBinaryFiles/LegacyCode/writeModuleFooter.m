function ok = writeModuleFooter(fid, footer)
    ok = true;
    try
        fwrite(fid, footer.length, 'int32');
        fwrite(fid, footer.identifier, 'int32');
        fwrite(fid, footer.binaryLength, 'int32');
        fwrite(fid, footer.binData, 'int8');
    catch
        ok = false;
    end