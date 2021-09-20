%==========================================================================
% CORAL AGE MODEL DEMO

% Helper script to demonstrate how to use the coral age modeling function
% psAgeModel.m using generic periodic data as the input for the algorithm. 

% The age model function is available at: https://github.com/lawmana/coralPSM.

% NOTE: The age model algorithm assumes there is cyclicity in the input
% data. For optimal use with coral geochemical data, the coral should be
% sampled at ~monthly resolution (12 samples/year). The annual cycle should
% be the dominant signal such that the age model algorithm can identify
% peaks/troughs in the raw data and interpolate to monthly-resolution.

% Author: A.E. Lawman, The University of Texas at Austin, 2020

% Citation: Lawman et al., (2020), Developing a coral proxy system model to
% compare, Paleoceanography and Paleoclimatology, doi: 10.1029/2019PA003836
%==========================================================================
%   AGE MODEL USER INPUTS:
%   dataIn      : the input data with a cyclic signal of known frequency
%   ppcIn       : points-per-cycle for the input data (i.e., the sampling resolution)
%   ppcOut      : points-per-cycle for the age modeled output 
%   peakLoc     : the index for each peak in each cycle (peakLoc must be between 1 and ppcOut)
%   troughLoc   : the index for each trough in each cycle (troughLoc must be between 1 and ppcOut)
%   numPeriods  : the estimated number of complete periods in the input data

%   AGE MODEL OUTPUTS:
%   dataOut     : the relative age modeled output [RelativeTime Data]
%   criticalPts : the indices for the local minima/maxima (troughs/peaks) in a cycle
%==========================================================================
close all; clc;

% dock figures
set(0,'DefaultFigureWindowStyle','docked')

%--------------------------------------
% GENERATE SAMPLE DATA TO AGE MODEL
%--------------------------------------
% sampling resolution for the input data (points-per-cyle)
% for sub-annually resolved coral data this is the estimated number of 
% samples per year (looking at the annual growth bands or the geochemistry)
ppcIn = 14;

% depth or sample number series with temporal uncertainty
noiseX = 0.5;
xIn = 1:1/ppcIn:20;
xnIn = xIn + (noiseX/ppcIn)*randn(size(xIn));

% sample sine wave data set with Gaussian noise (define amplitudeY to set 
% the amplitude of the sine wave)
noiseY = 0.1;
amplitudeY = 1.5;
dataIn = amplitudeY*sin(2*pi*xnIn);
dataIn = dataIn + noiseY*randn(size(dataIn));

%----------------------------------------------------------
%% EXAMPLE 1: no variable arguments
%----------------------------------------------------------
exampleType = 'Default settings (no variable arguments)';

% age model the input data
% criticalPts provides the indices for the local minima/maxima in dataIn
[dataOut, criticalPts] =  psAgeModel(dataIn);

% plot the raw data with the peaks/troughs identified by the age model
% plot the age modeled data and its corresponding climatology
figure(1); clf; hold on;
plotTimeSeries(exampleType, xIn, dataIn, criticalPts, dataOut)
text(0.025, 0.95, 'default ppcOut may not yield an an accurate climatology unless ppcIn = 12', 'units', 'normalized')

%----------------------------------------------------------
%% EXAMPLE 2: specify points-per-cycle (ppcIn and ppcOut)
%----------------------------------------------------------
exampleType = 'specify points-per-cycle (ppcIn and ppcOut)';

% Specify the number of ppc for the input data ('ppcIn') & the target
% temporal resolution for the age modeled output ('ppcOut'). If the annual 
% cycle is the dominant signal, set ppcOut = 12 to achieve monthly resolution
ppcOut = 12;

[dataOut, criticalPts] =  psAgeModel(dataIn, 'ppcIn', ppcIn, 'ppcOut', ppcOut);

figure(2); clf; hold on;
plotTimeSeries(exampleType, xIn, dataIn, criticalPts, dataOut)

%----------------------------------------------------------
%% EXAMPLE 3: peak/trough assignment
%----------------------------------------------------------
exampleType = 'specify ppcIn, ppcOut, and the peak/trough assignment';

% Specify the calendar months that correspond to the peaks and troughs 
% in the input data (e.g., the climatological warmest and coldest months)
%1 = January, 2 = February, ..., 11 = November, 12 = December.
peakLoc = 6;        
troughLoc = 11;      

[dataOut, criticalPts] =  psAgeModel(dataIn, 'ppcIn', ppcIn, 'ppcOut', ppcOut, 'peakLoc', peakLoc, 'troughLoc', troughLoc);

figure(3); clf; hold on;
plotTimeSeries(exampleType, xIn, dataIn, criticalPts, dataOut)

%----------------------------------------------------------
%% EXAMPLE 4: constrain the number of periods
% (e.g., the estimated number of years in the data)
%----------------------------------------------------------
exampleType = 'specify ppcIn, ppcOut, the peak/trough assignment, and estimated # of periods';

% Note: This constraint is often not necessary for data with clear
% cyclicity, but may be necessary for noisy data. 

% For corals, the number of cycles (years) constraint can be estimated by 
% counting the number of annual density bands that are visible in an x-ray 
% image of the cross-sectional slabs.

% estimated number of years in the input data
numYrs = 20;

[dataOut, criticalPts] =  psAgeModel(dataIn, 'ppcIn', ppcIn, 'ppcOut', ppcOut, 'peakLoc', peakLoc, 'troughLoc', troughLoc, 'numPeriods', numYrs);

figure(4); clf; hold on;
plotTimeSeries(exampleType, xIn, dataIn, criticalPts, dataOut)

%----------------------------------------------------------
%% EXAMPLE 5: constrain number of tie-points
% (e.g., assign just peaks)
%----------------------------------------------------------
exampleType = 'specify ppcIn, ppcOut, peak assignment, estimated # of periods, and # tie points';

% Note: This constraint is often not necessary for data with clear
% cyclicity, but may be necessary for noisy data. 

% For corals, the number of cycles (years) constraint can be estimated by 
% counting the number of annual density bands that are visible in an x-ray 
% image of the cross-sectional slabs.

% estimated number of years in the input data
numYrs = 20;

[dataOut, criticalPts] =  psAgeModel(dataIn, 'ppcIn', ppcIn, 'ppcOut', ppcOut, 'peakLoc', peakLoc, 'numTiePoints', 1, 'numPeriods', numYrs);

figure(5); clf; hold on;
plotTimeSeries(exampleType, xIn, dataIn, criticalPts, dataOut)


%% HELPER FUNCTIONS

% plotting
function [] = plotTimeSeries(exampleType, xIn, dataIn, criticalPts, dataOut)
    % climatology
    months = 1:12;
    clim = climatology(dataOut(:,2));
    [climMax, iMax] = max(clim);
    [climMin, iMin] = min(clim);
    
    % peaks/troughs
    peaks = criticalPts(dataIn(criticalPts) > 0);
    troughs = criticalPts(dataIn(criticalPts) <= 0);
    
    % plot formatting
    lw = 1.5;
    markerSz = 10;

    ax1 = subplot(2,4,1:3); hold on;
    title('Raw Data')
    plot(xIn, dataIn, '-k', 'linewidth', lw);
    plot(xIn(peaks), dataIn(peaks), 'ro', 'linewidth', lw, 'markersize', markerSz);
    plot(xIn(troughs), dataIn(troughs), 'bo', 'linewidth', lw, 'markersize', markerSz);
    xlabel('depth or sample #')
    legend('raw data', 'peaks', 'troughs', 'location', 'northwest', 'orientation', 'horizontal')
    
    ax2 = subplot(2,4,5:7); hold on;
    title(['Age Modeled Data: ' exampleType])
    plot(dataOut(:,1), dataOut(:,2), '-k', 'linewidth', lw)
    xlabel('relative year')

    ax3 = subplot(2,4,8); hold on;
    title('Climatology')
    plot(months, clim, '-ok', 'linewidth', lw, 'MarkerFaceColor', 'k')
    plot(iMax, climMax, 'r^', 'markersize', markerSz, 'MarkerFaceColor', 'r')
    plot(iMin, climMin, 'bv', 'markersize', markerSz, 'MarkerFaceColor', 'b')
    xlabel('calendar month')
    xlim([1 12])
    legend('climatology', 'avg. peak', 'avg. trough', 'location', 'southwest')
    
    % link axes for equal y-limits
    linkaxes([ax1 ax2 ax3],'y')
    
end

% climatology
function [clim] = climatology(data)
    
    clim = zeros(12,1);
    
    for month=1:12
        clim(month) = nanmean(data(month:12:end));
        
    end
    
end

