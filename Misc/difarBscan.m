function [ang freqs mag Output degBins freqBins OmEWx] = difarBscan(Om, EW, NS, fftLength, sampleRate, freqRange)
% [OutBRR degBins freqBins] = difarBscan(Om, EW, NS, sampleRate)
% Find the magnetic bearing of a sound source given the demultiplexed output
% of a difar sonobuoy, type AN-SSQ-53B or 53D
%
% The basis for this program is explained in a paper, "Relationship of
% Underwater Acoustic Intensity Measurements to Beamforming", by Gerald D'Spain,
% Canadian Acoustics, Proceedings Issue, September 1994, pp. 157-158
%
% This version plots the ambiguity surface for a Bartlett (linear) beamformer
% as sound frequency versus azimuth to the source
%
% read the three files, Omnidirectional pressure (Om), East-West Velocity (EW)
% and North-South Velocity (NS)
% as16 bit binary files with no header, 'Fso' sampling rate
%
% M. McDonald, 9/7/98
%
% Changes made by B. S. Miller, Australian Antarctic Division, July 2011
% FFTLength is the time in seconds for one spectrogram/pwelch PSD slice
% sampleRate is the sampling rate in Hz for the omni, ew, and ns signals
% freqRange is a 1x2 vector containing the min and max frequencies over
% which bearings will be computed.
%   -Cleaned up code 
%     -Renamed many variables and removed unused variables and commented
%      lines
%     -separated computation from plotting
%   -Changes to plotting
%     -Made plotting optional so that program can be run "batch" style
%     -Added compass plot showing bearings from all frequencies. [Removed
%        this because it was misleading as the bearings were not corrected
%        for magnetic deviation.
%     -Plot bearing-frequency surface, spectrogram, and a compass as 
%        subplots in the same figure window to reduce clutter on the
%        screen.
%   -Further efforts to improve "batch" style automation
%     -Function now takes as input vectors of omnidirectional pressure and
%       ns and ew velocity.
%     -Function now outputs the ambiguity surface as well as frequency and
%       angle bins
makePlots = 0;

if nargin < 0; freqRange = [-inf inf]; end
% convert fftLength in seconds, to fftLength in samples (make sure its an
% even number)
fftLength = fix(fftLength * sampleRate);
fftLength = min([fftLength length(Om)]);
if(2*fix(fftLength/2)~=fftLength)
    fftLength = fftLength - 1;
end

% The routines sppowr and spcros are from "Signal Processing Algorithms in Matlab"
% by S.D. Stearns and R. A. David, Prentice Hall 1996, but the Matlab signal
% processing routines (psd & csd) would probably work also
% NDFT=512;   % FFT length in samples

iwindo=4;   % Hanning window
ovrlap=0.5; % 50 percent overlap

% compute the autospectra of each component
OmS = pwelch(Om', hann(512), 256, 512);
EWS = pwelch(EW', hann(512), 256, 512);
NSS = pwelch(NS', hann(512), 256, 512);
% OmS = sppowr(Om', fftLength, iwindo, ovrlap);
% EWS = sppowr(EW', fftLength, iwindo, ovrlap);
% NSS = sppowr(NS', fftLength, iwindo, ovrlap);

% compute the cross spectra needed
OmEWx = cpsd(Om', EW', hann(512), 256, 512);
OmNSx = cpsd(Om', NS', hann(512), 256, 512);
EWNSx = cpsd(EW', NS', hann(512), 256, 512);
% compute the cross spectra needed
% OmEWx = spcros(Om', EW', fftLength, iwindo, ovrlap);
% OmNSx = spcros(Om', NS', fftLength, iwindo, ovrlap);
% EWNSx = spcros(EW', NS', fftLength, iwindo, ovrlap);

nfbins=(fftLength/2)+1;
fstep=(sampleRate/2)/nfbins;
freqBins=((1:nfbins)*fstep)-(fstep/2); % the center frequencies of each bin

%the frequency  bins below 15 Hz are discarded
fbinmin=ceil(15/fstep); 
fbinmin = max([fbinmin fix(freqRange(1)/fstep)]);

% the top 10 percent of the bins are discarded because they are in the
% rollof of the resampling filter
fbinmax=floor(nfbins-(0.1*nfbins));
fbinmax = min([fbinmax fix(freqRange(2)/fstep)]);

%the steering vectors are a 1X3 matrix for each step in azimuth (theta)
degstep=2;  % step in degrees
step=2*pi/180;  % step in radians
nazsteps=360/degstep;   % number of steps in azimuth
degBins = 1:degstep:360;
OutB=ones([nazsteps,nfbins]);

outputType = 'Bartlett'; % Can be 'Bartlett' or 'MVDR'
for az=1:1:nazsteps,
    theta = (az-1)*step;
    svec=[0.5, 0.5*sin(theta), 0.5*cos(theta)];
    for f = fbinmin:fbinmax, % for each frequency bin

        % form the 3X3 Data Cross Spectral Matrix, ignoring the conversion factors rho-c
        Sxm = [OmS(f),  OmEWx(f), OmNSx(f);...
              OmEWx(f), EWS(f),   EWNSx(f);...
              OmNSx(f), EWNSx(f), NSS(f) ];
        
        switch outputType  
            case 'Bartlett' %(McDonald version)
                OutB(az,f)=svec*Sxm*svec'; 
            case 'MVDR'     % MVDR beamformer (Thode modification)
                OutB(az,f)=1./(svec*inv(Sxm)*svec');
        end
     end
end

OutBR=real(OutB);  %real part of complex number

nazsteps2=nazsteps/2;

% Not so sure about this reversal, all my bearings are 180 degrees off...
%   BSM - 2012-January
% result is direction of energy flow, to get bearing shift the azimuths 180 degrees
% OutBRR = ones(nazsteps,nfbins);
% OutBRR(1:nazsteps2,:) = OutBR(nazsteps2+1:nazsteps,:);
% OutBRR(nazsteps2+1:nazsteps,:) = OutBR(1:nazsteps2,:);

freqhi=freqBins(fbinmin:fbinmax);

Output=log10(abs(OutBR(:,fbinmin:fbinmax)));

[mag,maxOutIx] = (max(Output)); % Maximum
ang = degBins(maxOutIx);
freqs = freqBins(fbinmin:fbinmax);

if makePlots
    iwindo=4;   % Hanning window
    ovrlap=0.5; % 50 percent overlap
    
    % find the max intensity azimuths for the frequencies with the most energy
    numFreqs = 4;
    [maxOut,maxOutIx] = (max(Output)); % Maximum
    [maxSorted sortIx] = sort(maxOut,'descend');
    angSorted = degBins(maxOutIx(sortIx(1:numFreqs)));
    freqSorted = freqhi(sortIx(1:numFreqs));

    [Y1,X1]=meshgrid(freqhi,degBins);
    
    
    hFig = figure;
    pos = get(hFig,'position');
    set(hFig,'position',[pos(1:2) 1000 600]);
    grid on;
    bearingFreq = subplot(1,2,2);
    pcolor(X1,Y1,Output);
    grid on;
    set (gca,'XTick', 45:45:315);
    shading flat;
%     vertLabel;
    xlabel('source bearing, magnetic');
    ylabel('frequency, Hz');

    hold on;
    colorbar;

    % figure(fh);
    plot(degBins(maxOutIx),freqhi,'xm');
    plot(angSorted,freqSorted,'o');

    %Below is JG addition to plot quick spectrogram of sound clip
    timeFreq = subplot(1,2,1);
    windo = min([sampleRate length(Om)]);

    spectrogram(Om, windo ,fix(windo*.9),windo,sampleRate,'yaxis');
    ylim([freqhi(1) freqhi(end)]);
    [cmin,cmax]=caxis;
    caxis([cmin/1.5,cmax]); %JG addition--to restrict colormap so values lower than cmin/1.5 are all one color...makes more intense sounds standout.
    colorbar;

    linkaxes([bearingFreq, timeFreq],'y');
    
    % The code below has been commented in order to reduce confusion. It
    % does not account for magnetic deviation, since this is done at a
    % later stage of processing.
    % A polar plot with arrows pointing in the direction of the sound. Arrow
    % length is proportional to the unlogged energy.
%     subplot(1,3,3);
    % polar(degBins(maxOutIx)*pi/180,maxOut);
%     [x y] = pol2cart(degBins(maxOutIx)*pi/180,10.^maxOut);
%     compass(x,y);
    % polar(repmat(()*pi/180,size(Output,2),1)',Output.*(Output>(maxOut - 0.1)))
end
