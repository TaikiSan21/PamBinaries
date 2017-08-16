function ok = writeClick(file, click)
% first work out the total length of the click. 
clickLen = getClickDataLen(click);
ok = true;
try
    fwrite(file, clickLen+8, 'int32');
    fwrite(file, 1000, 'int32');
    fwrite(file, click.timeMillis, 'int64');
    fwrite(file, clickLen-4, 'int32');
    fwrite(file, click.startSample, 'int64');
    fwrite(file, click.channelMap, 'int32');
    fwrite(file, click.triggerMap, 'int32');
    fwrite(file, click.type, 'int16');
    fwrite(file, click.flags, 'int32');
    fwrite(file, length(click.delays), 'int16');
    fwrite(file, click.delays, 'float');
    fwrite(file, length(click.angles), 'int16');
    fwrite(file, click.angles, 'float');
    fwrite(file, length(click.angleerrors), 'int16');
    fwrite(file, click.angleerrors, 'float');
    fwrite(file, click.duration, 'int16');
    maxval = max(abs(click.wave));
    fwrite(file, maxval, 'float');
    fwrite(file, click.wave*127/maxval, 'int8');
%     disp('click written ok')
catch
    ok = false
end