%==========================================================================
% CORAL ANALYTICAL ERROR EXAMPLE: VANUATU FOSSIL CORAL Sr/Ca 

% Author: A.E. Lawman, The University of Texas at Austin, 2020

% Generate realizations of fossil coral Sr/Ca data perturbed with Gaussian
% analytical error (instrument error).

% This examples uses Sr/Ca data from fossil coral 11-TM-S5 collected from 
% Tasmaloum, Vanuatu (lat/lon: 15.6S, 166.9E).
%--------------------------------------------------------------------------
% ANALYTICAL_ERROR USER INPUTS 
% see function documentation for more details
%  X               : the input data, here coral Sr/Ca data vs. depth in a
%  sigma           : the standard deviation for the analytical error term
%                    (assumed to be Gaussian). Multiply by 2 for ±2sigma if
%                    desired
%  numRealizations : the desired number of realizations

% ANALYTICAL_ERROR OUTPUTS:
%  Xn              : a length(X) x numRealizations matrix where each column
%                    is a realization of the input data perturbed with 
%                    Gaussian noise
%--------------------------------------------------------------------------
% ANALYTICAL_ERROR ALGORITHM CITATION:
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
%% USER INPUTS
%=================================
sigma = 0.012; % units (mmol/mol); ±2sigma
numRealizations = 10000;

%=================================
%% LOAD CORAL DATA
%=================================
% approximately monthly-resolved Sr/Ca data from fossil coral 11-TM-S5
% collected from Tasmaloum, Vanuatu (lat/lon: 15.6S, 166.9E)
% depth (in mm from the top of the core) vs. Sr/Ca
load 'Vanuatu_FossilCoral_11-TM-S5_SrCa.mat'

% fill missing data in the raw Sr/Ca series using linear interpolation (for
% illustrative plotting purposes to avoid issues with NaNs when using 'patch')
SrCa = fillmissing(SrCa,'linear');

%=================================
%% GENERATE Sr/Ca REALIZATIONS
%=================================
% generate realizations perturbed with Gaussian noise
realizations = analytical_error(SrCa, sigma, numRealizations);

% calculate the minimum and maximum values across all realizations at each
% depth (for plotting the error cloud), take min/max across the second
% dimension of the realization matrix
minSrCa = min(realizations, [], 2);
maxSrCa = max(realizations, [], 2);

%=================================
%% FIGURE GENERATION
%=================================
alpha = 0.5; % transparency for the patch
plotColor = [0 0 0.75]; % RGB value

% plot the original Sr/Ca input data with shaded cloud for the
% realizations perturbed with analytical_error
figure(1); clf; hold on;
title({[coralType ' Coral Sr/Ca ' coralID ' Perturbed with Analytical Error'] ...
    ['n = ' num2str(numRealizations) ' realizations']})

% realizations plotted using patch
patch([depth; flip(depth)], [minSrCa; flip(maxSrCa)], plotColor,'FaceAlpha',alpha, 'EdgeColor', 'none');

% original input data
plot(depth, SrCa, '.-', 'color', plotColor)

xlabel('depth from top of core (cm)')
ylabel('Sr/Ca (mmol/mol)')
legend('±2\sigma', 'Sr/Ca input', 'location', 'northwest')

% axis formatting
ax = gca;
ax.FontWeight = 'bold';
ax.XDir = 'reverse';
ax.YDir = 'reverse';
ax.XMinorTick = 'on';
