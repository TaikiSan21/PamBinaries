function click = createClick(timeMillis, sample, duration, wave);
click.timeMillis = timeMillis;
click.startSample = sample;
click.duration = duration;
click.wave = wave;
sz = size(wave);
click.nChan = min(sz);
chanMap = 0;
for i = 1:click.nChan
    chanMap = chanMap + 2^(i-1);
end
click.channelMap = chanMap;
click.triggerMap = chanMap;
click.type = 0;
click.flags = 0;
click.delays = [];
click.angles = [];
click.angleerrors = [];