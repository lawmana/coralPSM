function[Xn] = analytical_error(X, sigma, numRealizations)
% ANALYTICAL_ERROR(X, sigma, numRealizations)
%     ANALYTICALERROR adds Gaussian white noise to simulate
%     analytical / measurement error.
%        
%     Input Arguments:
%           1. X, vector (1D)
%           2. sigma (assumed precision of measurement/instrumental error)
%           3. nRealizations (number of noise realizations)
%     Output:
%           Returns a matrix (size length(X) x numRealizations), Xn includes 
%           analytic error due to measurement precision.
%          
%           Example of use:  generate 100 realizations of speleothem proxy
%           data with analytical error.
%           proxy_record = analytical_error(d18O, 0.1, 100)

%==========================================================================
    % input checks
    % X is a column vector
    if ~isvector(X)
        error('Incorrect user input: X must be a column vector')
    end
    
    % sigma and nsamples are scalars
    if ~isscalar(sigma)
        error('Incorrect user input: sigma must be a scalar')
    end
    
    if ~isscalar(numRealizations)
        error('Incorret user input: nsamples must be a scalar')
    end
    
    % ensure that X is a column vector
    if isrow(X)
        X = X';
    end
%==========================================================================
    
% generate vector/matrix of Gaussian noise and add to the input vector
    noise = sigma.*randn([size(X,1),numRealizations]);
    Xn    = repmat(X,1,numRealizations) + noise;
end