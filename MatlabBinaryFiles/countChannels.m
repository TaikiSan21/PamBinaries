function nC = countChannels(channelMap)
% function nC = countChannels(channelMap)
% count the numebr of set bits in the channel map
nC = 0;
j = 1;
for i = 1:32
    if (bitand(channelMap, j))
        nC = nC + 1;
    end
    j = j * 2;
end