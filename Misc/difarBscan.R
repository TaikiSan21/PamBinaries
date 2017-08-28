# Find the magnetic bearing of a sound source given the demultiplexed
# output of a difar sonobuoy, type AN-SSQ-53B or 53D

# The basis for this program is explained in a papaer, "Relationship of
# Underwater Aocustic Intensity Measurements to Beamforming" by Gerald
# D'Spain, Canadian Acoustics, Proceedings Issue, Sept. 1994, pp.157-158

# This version plots the ambiguity surface for a Bartlett (linear) 
# beamformer as sound frequency versus azimuth to the source

# read the three files, Omnidirectional pressure (Om), East-West Velocity (EW)
# and North-South Velocity (NS)
# as 16-bit binary files with no header, 'Fso' sampling rate

# M. McDonald, 9/7/98

# Changes made by B. S. Miller, Australian Antarctic Division, July 2011
# FFTLength is the time in seconds for one spectrogram/pwelch PSD slice
# sampleRate is the sampling rate in Hz for the omni, ew, and ns signals
# freqRange is a 1x2 vector containing the min and max frequencies over
# which bearings will be computed.
#   -Cleaned up code
#       -Renamed many variables and removed unused variables and commented lines
#       - Separated computation from plotting
#   -Changes to plotting
#       -Made plotting optional so that program can be run 'batch' style
#       -Added compass plot showing bearings from all frequencies. [Removed
#        this because it was misleading as the bearings were not corrected
#        for magnetic deviation.
#       -Plot bearing-frequency surface, spectrogram, and a compass as
#        subplots in the same figure window to reduce clutter on the screen.
#   -Further efforts to improve 'batch' style automatino
#       -Function now takes as input vectors of omnidirectional pressure and ns
#        and ew velocity
#       -Function now outputs the ambiguity surface as well as frequency and
#        angle bins

# This code adapted from Matlab to R by Taiki Sakai

difarBscan <- function(Om, EW, NS, fftLength, sampleRate, freqRange=c(-Inf, Inf)) {
    makePlots = FALSE
    # Convert fftLength in seconds, to fftLength in samples
    # Don't think this is needed for R version of PSD.
    
    fftLength <- trunc(fftLength * sampleRate)
    fftLength <- min(fftLength, length(Om))
    # # Make sure it is an even number
    fftLength <- fftLength - (fftLength %% 2)
    
    # The routines sppowr and spcros from "Signal Processing Algorithms in Matlab"
    # by S.D. Stearns and R. A. David, Prentice Hall 1996 were used in Matlab version.
    # Here using the sapa package.
    
    # Computes entire cross spectra matrix. S11, S12, S13, S22, S23, S33. Uses Hanning
    # window with overlap of .5.
    allSpectra <- SDF(matrix(c(Om, EW, NS), ncol=3), method='wosa', sampling.interval=1/sampleRate)
    
    OmS <- allSpectra[,1]
    EWS <- allSpectra[,4]
    NSS <- allSpectra[,6]
    
    OmEWx <- allSpectra[,2]
    OmNSx <- allSpectra[,3]
    EWNSx <- allSpectra[,5]
    
    nfbins <- (fftLength/2)+1
    fstep <- (sampleRate/2)/nfbins
    freqBins <- ((1:nfbins)*fstep)-(fstep/2) # the center of frequencies of each bin
    
    # The frequency bins below 15Hz are discarded
    fbinmin <- ceiling(15/fstep)
    fbinmin <- max(fbinmin, trunc(freqRange[1]/fstep))
    
    # The top 10 percent of the bins are discarded because they are in the
    # rolloff of the resampling filter
    fbinmax <- floor(nfbins-(.1*nfbins))
    fbinmax <- min(fbinmax, trunc(freqRange[2]/fstep))
    
    # The steering vectors are a 1x3 matrix for each step in azimuth (theta)
    degstep <- 2 # step in degrees
    step  <- degstep*pi/180 # step in radians
    nazsteps  <- 360/degstep
    degBins <- seq(1, 360, degstep)
    OutB <- matrix(1, nrow=nazsteps, ncol=nfbins)
    
    outputType <- 'Bartlett' # Can be 'Bartlett' or 'MVDR'
    for(az in 1:nazsteps) {
        theta <- (az-1)*step
        svec <- c(.5, .5*sin(theta), .5*cos(theta))
        for(f in fbinmin:fbinmax) { # for each frequency bin
            
            # form the 3x3 Data Cross Spectral Matrix, ignoring the conversion factors rho-c
            Sxm <- matrix(c(OmS[f], OmEWx[f], OmNSx[f],
                         OmEWx[f], EWS[f], EWNSx[f],
                         OmNSx[f], EWNSx[f], NSS[f]),
                         nrow=3, ncol=3)
            
            switch(outputType,
                'Bartlett' = {
                    OutB[az, f] <- svec%*%Sxm%*%svec
                },
                'MVDR' = {
                    OutB[az, f] <- 1/(svec%*%solve(Sxm)%*%svec)
                })
        }
    }
    OutBR <- Re(OutB)
    
    nazsteps2 <- nazsteps/2
    
    # Not so sure about this reversal, all my bearings are 180 degrees off...
    #   BSM - 2012-January
    # result is direction fo energy flow, to get bearing shift the azimuths 180 deg
    # OutBRR <- matrix(1, nrow=nazsteps, ncol=nfbins)
    # OutBRR[1:nazsteps2,] <- OutBR[(nazsteps2+1):nazsteps,]
    # OutBRR[(nazsteps2+1):nazsteps,] <- OutBR[1:nazsteps2,]
    
    freqhi <- freqBins[fbinmin:fbinmax]
    
    Output <- log10(abs(OutBR[,fbinmin:fbinmax]))
    
    # Max of each column, and index of that max
    mag <- apply(Output, 2, max)
    maxOutIx <- apply(Output, 2, which.max)
    ang <- degBins[maxOutIx]
    # This seems weird?
    freqs <- freqBins[fbinmin:fbinmax]
    return(list(ang=ang, freqs=freqs, mag=mag, Output=Output,
                degBins=degBins, freqBins=freqBins))
}
    