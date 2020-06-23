%==========================================================================
% CORAL AGE MODEL DEMO

% Script to demonstrate how to use the coral age modeling function
% psAgeModel.m available at: https://github.com/lawmana/coralPSM

% NOTE: The age model algorithm assumes there is cyclicity in the input
% data. For corals, it is assumed that the coral was sampled at
% approximately monthly resolution (12 points per year), and that the
% annual cycle is the dominant signal.

% A.E. Lawman
% The University of Texas at Austin, 2020

% Citation: Lawman et al., (2020), Developing a coral proxy system model to
% compare, Paleoceanography and Paleoclimatology, doi: <enter here>
%==========================================================================
clear all; close all; clc;

% dock figures
set(0,'DefaultFigureWindowStyle','docked')

%--------------------------------------
% GENERATE SAMPLE DATA
%--------------------------------------
% sampling resolution for the input data (points-per-cyle)
ppcIn = 12;

% depth or sample number series with temporal uncertainty
noiseX = 0.3;
xIn = 1:1/ppcIn:20;
xnIn = xIn + (noiseX/ppcIn)*randn(size(xIn));

% sample data set with Gaussian noise
noiseY = 0.1;
dataIn = sin(2*pi*xnIn);
dataIn = dataIn + noiseY*randn(size(dataIn));

%----------------------------------------------------------
%% EXAMPLE 1: no variable arguments
%----------------------------------------------------------
% age model the input data
[dataOut, criticalPts] =  psAgeModel(dataIn);

% plot the original data with the peaks/troughs identified by the age model
% plot the age modeled data and its corresponding climatology
figure(1); clf; hold on;
plotTimeSeries(xIn, dataIn, criticalPts, dataOut)

%----------------------------------------------------------
%% EXAMPLE 2: points-per-cycle (ppc)
%----------------------------------------------------------
% Specify the number of ppc for the input data ('ppcIn') & the target
% temporal resolution for the age modeled output ('ppcOut')
ppcOut = 12;

[dataOut, criticalPts] =  psAgeModel(dataIn, 'ppcIn', ppcIn, 'ppcOut', ppcOut);

figure(2); clf; hold on;
plotTimeSeries(xIn, dataIn, criticalPts, dataOut)

%----------------------------------------------------------
%% EXAMPLE 3: assign calendar months to the peaks/troughs 
%----------------------------------------------------------
% Specify the calendar months that correspond to the peaks and troughs 
% in the input data (e.g., the climatological warmest and coldest months)
%1 = January, 2 = February, 11 = November, 12 = December.
peakLoc = 2;        % February
troughLoc = 8;      % August

[dataOut, criticalPts] =  psAgeModel(dataIn, 'ppcIn', ppcIn, 'ppcOut', ppcOut, 'peakLoc', peakLoc, 'troughLoc', troughLoc);

figure(3); clf; hold on;
plotTimeSeries(xIn, dataIn, criticalPts, dataOut)

%----------------------------------------------------------
%% EXAMPLE 4: constrain the number of periods
% (e.g., the estimated number of years in the data)
%----------------------------------------------------------
% Note: This constraint is often not necessary for data with a clear annual
% cycle, but may be necessary for noisy data. For corals, the number of
% cycles (years) constraint can be estimated by counting the number of
% annual density bands that are visible in an x-ray image.

% estimated number of years in the input data
numYrs = 20;

[dataOut, criticalPts] =  psAgeModel(dataIn, 'ppcIn', ppcIn, 'ppcOut', ppcOut, 'peakLoc', peakLoc, 'troughLoc', troughLoc, 'numPeriods', numYrs);

figure(4); clf; hold on;
plotTimeSeries(xIn, dataIn, criticalPts, dataOut)

%% HELPER FUNCTIONS

function [] = plotTimeSeries(xIn, dataIn, criticalPts, dataOut)
    months = 1:12;
    clim = climatology(dataOut(:,2));
    
    lw = 1.5;

    subplot(2,4,1:3); hold on;
    title('Input Data & Critical Points')
    plot(xIn, dataIn, '-k', 'linewidth', lw);
    plot(xIn(criticalPts), dataIn(criticalPts), 'ro', 'linewidth', lw);

    subplot(2,4,5:7); hold on;
    title('Age Modeled Data')
    plot(dataOut(:,1), dataOut(:,2), '-k', 'linewidth', lw)

    subplot(2,4, 8); hold on;
    title('Climatology')
    plot(months, clim, '-ok', 'linewidth', lw, 'MarkerFaceColor', 'k')
    
end

function [clim] = climatology(data)
    
    clim = zeros(12,1);
    
    for month=1:12
        clim(month) = nanmean(data(month:12:end));
        
    end
    
end

