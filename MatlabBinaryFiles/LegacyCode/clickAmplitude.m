function amp = clickAmplitude(clickWave, hSens, gain, adcPeakPeak)
sz = size(clickWave);
nChan = sz(2);
amp = zeros(1,nChan);
for i = 1:nChan
    mV = max(abs(clickWave(:,i)));
    dbV = 20*log10( mV/2*adcPeakPeak );
    amp(i) = dbV - (hSens + gain);
end