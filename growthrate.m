function[perturbedData, depthVec] = growthrate(X, Y, lagcoefs, growthRateSigma)
% GROWTHRATE(X,Y,lagcoefs, growthRateSigma, growthRateMu)
% GROWTHRATE perturbs the independent vector of a data set using an
% autoregressive (AR) model. 
% INPUT:
    % X is a 1-D vector (independent variable)
    % Y is a 1-D vector (dependent variable)
    % LAGCOEFS is a 1-D vector of coefficients for the autoregressive model
    %   E.x. [lag1coef lag2coef, lag3coef ... lagNcoef]
    % GROWTHRATESIGMA is a scaler input for the standard deviation of the
    %   AR model
% OUTPUT:
    % PERTURBEDDATA is 1-D vector with the perturbed data following a
    % specified autoregressive model based on user lag coefficients, mean,
    % and standard deviation
    % DEPTHVEC is a 1-D vector
    
%==========================================================================    
    % input checks
    if ~isscalar(growthRateSigma)
        error('Incorrect user input: growthRateSigma must be a scalar')
    end

    if ~isvector(X)
        error('Incorrect user input: X must be a vector')
    end

    if ~isvector(Y)
        error('Incorrect user input: Y must be a vector')
    end

    if ~isvector(lagcoefs)
       error('Incorrect user input: lagcoefs must be a vector')
    end

    % ensure that X, Y and lagcoefs are column vectors
    if isrow(X)
       X = X';
    end

    if isrow(Y)
       Y = Y';
    end

    if isrow(lagcoefs)
       lagcoefs = lagcoefs';
    end
    
    % ensure that X and Y are the same size
    if ~(length(X) == length(Y))
        error('Incorrect user input: X and Y must be the same length');
    end
%========================================================================== 

    if ~exist('divPerRealization','var')
        divPerRealization = 12;
    end
    
    meanXWidth = mean(diff(X));
    
    n = length(X);
    numYears = ceil(n/divPerRealization);

    % AR model generated using a set of given lag coefficients (lagcoefs)
    Mdl = arima('Constant',0,'AR',num2cell(lagcoefs),'Variance',growthRateSigma*growthRateSigma);

    % evaluate the AR model and make the mean 1
    [perturbations,~] = simulate(Mdl, numYears);
    perturbations = perturbations + (1 - mean(perturbations));

    % interpolate the AR model to the same # data points as the
    % original depth series 
    perturbations_interp = interp1(1:numYears, perturbations, linspace(1, n/divPerRealization, n), 'linear'); % double check n/div vs numYears

    % adjust the Y-interp vector such that the average sampling resolution is
    % equal to the original vector
    perturbations_interp = meanXWidth * perturbations_interp/mean(perturbations_interp); % change to not use meanXWidth

    % generate monotonically increasing depth vector 
    depthVec = zeros(1,n);
    depthVec(1) = perturbations_interp(1);

    for j = 2:n
        depthVec(j) = depthVec(j-1) + perturbations_interp(j);
    end

    % interpolate the original data to the new perturbed
    % (compressed and stretched) depth vector
    perturbedData = interp1(X,Y, depthVec, 'linear');
    perturbedData = perturbedData';
end