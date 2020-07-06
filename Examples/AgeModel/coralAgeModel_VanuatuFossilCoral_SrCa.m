%==========================================================================
% CORAL AGE MODELING EXAMPLE: VANUATU FOSSIL CORAL Sr/Ca 

% Author: A.E. Lawman, The University of Texas at Austin, 2020

% Generate the relative age model for fossil coral 11-TM-S5 collected from
% Tasmaloum, Vanuatu (lat/lon: 15.6S, 166.9E).
%--------------------------------------------------------------------------
% AGE MODEL USER INPUTS:
%   dataIn      : the input data with a cyclic signal of known frequency; 
%                 here it is the raw coral Sr/Ca data in the depth domain
%                 with the annual cycle emerging as the dominant cyclic
%                 signal of known frequency. 
%   ppcIn       : points-per-cycle for the input data (i.e., the sampling resolution);
%                 here ppcIn represents the estimated # of coral samples per annual growth band
%   ppcOut      : points-per-cycle for the age modeled output
%                 set ppcOut = 12 for monthly-resolved ouput
%   peakLoc     : the index for each peak in each cycle (peakLoc must be between 1 and ppcOut)
%                 here this is the climatological warmest month
%   troughLoc   : the index for each trough in each cycle (troughLoc must be between 1 and ppcOut)
%                 here this is the climatological coolest month
%   numPeriods  : the estimated number of complete periods in the input data
%                 here this is the estimated number of years of coral
%                 geochemical data (based on the # of growth bands or # of
%                 annual cycles that appear in the raw geochemical data)

% AGE MODEL OUTPUTS:
%   dataOut     : the relative age modeled output [RelativeTime Data] stored as a Nx2 matrix
%   criticalPts : the indices for the local minima/maxima (troughs/peaks) in a cycle
%--------------------------------------------------------------------------
% AGE MODEL ALGORITHM CITATION:
% Lawman, A.E., Partin, J.W., Dee, S.G., Casadio, C.A., Di Nezio, P., 
% & Quinn, T.M. (2020). Developing a coral proxy system model to compare 
% coral and climate model estimates of changes in paleo-ENSO variability. 
% Paleoceanography and Paleoclimatology, 35, e2019PA003836. 
% doi: 10.1029/2019PA003836.

% CORAL DATA CITATION:
% Lawman, A.E., Quinn, T.M., Partin, J.W., Thirumalai, K., Taylor, F.W., 
% Wu, C.C., Yu, T.L., Gorman, M.K. & Shen, C.C. 2020. A Century of 
% Reduced ENSO Variability During the Medieval Climate Anomaly. 
% Paleoceanography and Paleoclimatology, 35, e2019PA003742. 
% doi: 10.1029/2019PA003742.

% NOAA Archive Link: https://www.ncdc.noaa.gov/paleo-search/study/28590
%==========================================================================
close all; clc;

%=================================
%% USER INPUTS FOR AGE MODELING
%=================================
%---Resolution inputs---
ppyIn = 12;         % estimated number of samples per year (points-per-year)
ppyOut = 12;        % temporal resolution (points-per-year) for the age-modeled data

%---Calendar month assignment for geochemical peaks/troughs---
% based on instrumental observations (e.g., HadISST; Rayner et al., [2003])
warmestMonth = 2;   % climatological warmest month at Vanuatu (February)
coolestMonth = 8;   % climatological coolest month at Vanuatu (August)

%---Estimated number of years of data---
% based on the number of annual growth bands in the x-ray image

% NOTE: this parameter is often not needed for data with a clear annual
% cycle, but is provided here for illustrative purposes
numYrs = 104;

%=================================
%% LOAD CORAL DATA
%=================================
% approximately monthly-resolved Sr/Ca data from fossil coral 11-TM-S5
% collected from Tasmaloum, Vanuatu (lat/lon: 15.6S, 166.9E)
% depth (in mm from the top of the core) vs. Sr/Ca
load 'Vanuatu_FossilCoral_11-TM-S5_SrCa.mat'

%=================================
%% GENERATE RELATIVE AGE MODEL
%=================================

% Step 1: invert the Sr/Ca data to reflect the inverse relationship between 
% temperature and coral Sr/Ca
% this will ensure that the identified peaks in the geochemical data with
% be Sr/Ca minima, and the identified trougs will be Sr/Ca maxima
SrCa_flip = -1 .* SrCa;

% Step 2: call psAgeModel.m including the user inputs defined above and the
% inverted Sr/Ca data
[ageMdl, criticalPts] = psAgeModel(SrCa_flip, 'ppcIn', ppyIn, 'ppcOut', ...
    ppyOut,'peakLoc', warmestMonth, 'troughLoc', coolestMonth, ...
    'numPeriods', numYrs);

% Step 3: invert the age modeled Sr/Ca output
ageMdl(:,2) = -1 .* ageMdl(:,2);

% calculate climatology (average seasonal cycle) for the monthly resolved
% age modeled time series
clim = climatology(ageMdl(:,2));

% calculate # of years of age modeled geochemical data
yrs = round(length(ageMdl)/ppyOut, 2);

%=================================
%% FIGURE GENERATION
%=================================
% fill missing data in the raw Sr/Ca series using linear interpolation (for
% illustrative plotting purposes only)
SrCa = fillmissing(SrCa,'linear');

figure(1); clf; hold on;
%---raw Sr/Ca data---
ax1= subplot(2, 3, 1:2); hold on;
title(['Raw ' coralType ' Coral Sr/Ca: ' coralID])
plot(depth, SrCa, '.-k')
plot(depth(criticalPts), SrCa(criticalPts), 'ro')
xlabel('depth from top of core (cm)')
ylabel('Sr/Ca (mmol/mol)')
legend('data', 'annual peaks & troughs', 'location', 'northwest')
ax1.XLim = [0 depth(1)];
ax1.XDir = 'reverse';
ax1.YDir = 'reverse';
ax1.XMinorTick = 'on';

%---text summary of the age model user inputs---
ax2 = subplot(2, 3, 3); hold on;
text(-0.2, 0.65, {'Summary of user inputs:', ...
    ['1. Estimated sampling resolution (samples/yr) for the raw data: ' num2str(ppyIn)], ...
    ['2. Temporal resolution for the age model (samples/yr): ' num2str(ppyOut)], ...
    ['3. Sr/Ca minima assignment (warmest month): ' num2str(warmestMonth)], ...
    ['4. Sr/Ca maxima assignment (coolest month): ' num2str(coolestMonth)], ...
    ['5. Estimated # of growth bands (years): ' num2str(numYrs)]}, 'units', 'normalized')
axis off

%---age modeled Sr/Ca data---
ax3= subplot(2, 3, 4:5); hold on;
title('Monthly Coral Sr/Ca (Age Modeled)')
plot(ageMdl(:,1), ageMdl(:,2), '.-k')
xlabel('relative year')
ylabel('Sr/Ca (mmol/mol)')
text(0.025, 0.95, ['# years = ' num2str(yrs)], 'units', 'normalized')
ax3.XLim = [0 ageMdl(end,1)];
ax3.YDir = 'reverse';
ax3.XMinorTick = 'on';

%---age modeled climatology---
ax4 = subplot(2, 3, 6); hold on;
title('Climatology (Age Modeled)')
plot(1:12, clim, 'o-k', 'MarkerFaceColor', 'k')
plot(warmestMonth, clim(warmestMonth), 'ro', 'MarkerSize', 10)
plot(coolestMonth, clim(coolestMonth), 'ro', 'MarkerSize', 10)
xlabel('calendar month')
ylabel('Sr/Ca (mmol/mol)')
ax4.XLim = [1 12];
ax4.YDir = 'reverse';
text(0.025, 0.90, {'Month assignment for peaks/troughs:', ['Sr/Ca minima: ' num2str(warmestMonth)], ...
    ['Sr/Ca maxima: ' num2str(coolestMonth)]}, 'units', 'normalized')

linkaxes([ax1 ax3 ax4], 'y')
%=================================
%% HELPER FUNCTIONS
%=================================
% climatology
function [clim] = climatology(data)
    
    clim = zeros(12,1);
    
    for month=1:12
        clim(month) = nanmean(data(month:12:end));
        
    end
    
end
