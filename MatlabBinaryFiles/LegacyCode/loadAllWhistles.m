function wsls = loadAllWhistles(bdFolder, fRange, fs, fftLen)
if nargin == 0
    bdFolder = 'C:\WhistleClassifier\BinaryBySepecies'
end
binaryEnds = {'pgdf', 'pgdx'}

wsls = [];

mask = ['\WhistlesMoans*.' binaryEnds{1}]
binaryFiles = dirsub(bdFolder, mask);
for f = 1:length(binaryFiles)
    aFile = binaryFiles(f).name;
    [tones, fileHead, fileFoot, modHead, modFoot] = loadWhistleFile(aFile);
    n1 = numel(tones);
    if nargin == 4
        tones = selectRange(tones, fRange, fs, fftLen);
    end
    disp(sprintf('%d of %d Whistle contours loaded from %s', ...
        length(tones), n1, aFile));
    wsls = [wsls tones];
end
%     eval(sprintf('save contours_%s.mat wsls', subFolders(i).name));
%     break
end

function selTones = selectRange(tones, fRange, fs, fftLen)
if isempty (tones)
    selTones = [];
    return;
end
  want = zeros(1, numel(tones));
  for i = 1:numel(tones)
      want(i) = wantTone(tones(i), fRange, fs, fftLen);
  end
  selTones = tones(find(want));
end

function want = wantTone(tone, fRange, fs, fftLen)
f = tone.contour*fs/fftLen;
fr = minmax(f);
if (fr(2) < fRange(1))
    want = 0;
elseif (fr(2) > fRange(2))
    want = 0;
else
    want = 1;
end
end