function [ts, criticalPts] = psAgeModel(y,varargin)
%psAgeModel Periodic signal tuned age model.
%   This function assumes that the largest signal in the given data is of
%   constant frequency and is distorted in time. The data will be
%   interpolated to even sampling in time by using the peaks and troughs of
%   the periodic signal as reference.

    %% parse input
    p = inputParser;

    addRequired(p, 'y', @(x) isnumeric(x) && isvector(x));
    addParameter(p, 'ppcIn', 0, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'ppcOut', 0, @isnumeric);
    addParameter(p, 'numPeriods', 0, @isnumeric);
    addParameter(p, 'peakLoc', 0, @isnumeric);
    addParameter(p, 'troughLoc', 0, @isnumeric);

    parse(p,y,varargin{:});

    ppcIn = p.Results.ppcIn;
    ppcOut = p.Results.ppcOut;
    peakLoc = p.Results.peakLoc;
    troughLoc = p.Results.troughLoc;
    numPeriods = p.Results.numPeriods;
    
    %% Initialize input
    
    % remove NaNs
    y = y(~isnan(y));
    
    % set y to a column vector
    if isrow(y)
        y = y';
    end
    
    % check length of y
    if length(y) < 10
        error('input data is too short in length');
    end
    
    %% Points per cycle (input): ppcIn
    
    % if no value is given for ppcIn then use a periodogram to find an approximation.
    if ppcIn == 0
        [pxx, w] = periodogram(y);
        [~,i] = max(pxx);
        ppcIn = 2*pi/w(i);
    end

    if ppcIn < 6
        error('ppc_in (%f) is too small',ppcIn);
    end

    %% Points per cycle (output): ppcOut
    
    % if no value is given for ppcOut then default to ppcIn rounded down to
    % the nearest integer
    if ppcOut == 0
        ppcOut = floor(ppcIn);
    end

    % error if output points per cycle is greater than input points per
    % cycle
    if ppcOut > ppcIn
        warning('input sampling rate (%f) < output sampling rate (%f).',ppcIn,ppcOut);
    end
    
    %% Peak location
    
    % set default peak location to maximum of seasonal cycle
    if peakLoc == 0
        [~,peakLoc] = max(periodicMean(y,round(ppcIn)));
    end

    if peakLoc > ceil(ppcIn)
        error('peak location (%f) is larger than points per cycle (%f).',peakLoc,ppcIn)
    end
    
    if peakLoc < 1
        error('peak location must be greater than or equal to 1. It was %f.',peakLoc);
    end
    
    %% Trough location
    
    % set default trough location to minimum of seasonal cycle
    if troughLoc == 0
        [~,troughLoc] = min(periodicMean(y,round(ppcIn)));
    end

    if troughLoc > ceil(ppcIn)
        error('trough location (%f) is larger than points per cycle (%f).',troughLoc,ppcIn)
    end
    
    if troughLoc < 1
        error('trough location must be greater than or equal to 1. It was %f.',troughLoc);
    end

    %% Peak to trough (p2t) and trough to peak (t2p) spacing
    
    if troughLoc > peakLoc
        p2t = troughLoc - peakLoc;
        t2p = ppcOut - p2t;
    else
        t2p = peakLoc - troughLoc;
        p2t = ppcOut - t2p;
    end
    
    %% Target number of periods (cycles)
    
    if numPeriods < 0
        error('number of years cannot be negative');
    elseif numPeriods == 0
        % default number of cycles based on the length of the input data and
        % ppcIn
        numPeriods = length(y)/ppcIn;
    end
    
    %% Lowpass filter
    
    % store initial input vector & apply low pass filter
    yInit = y;
    y = lowpass(y, 2, ppcIn);
       
    
    %% Find peaks
    
    % peak location (l) and prominence (p)
    [~,pks.l,~,pks.p] = findpeaks(y);

    currentLoc = -ppcIn/2;
    selectedPeaks = zeros(size(y));
    iterations = 0;
    runningValue = 0;
    while true
        iterations = iterations + 1;
        [currentLoc,value] = chooseNextPeak(currentLoc,pks.l,pks.p,ppcIn);
        selectedPeaks(currentLoc) = true;
        runningValue = runningValue + value;
        
        if isempty(currentLoc)
            break;
        end
        
        if currentLoc > (length(y) - 1.2*ppcIn)
            [currentLoc,value] = chooseNextPeak(currentLoc,pks.l,pks.p,ppcIn);
            
            if ~isempty(currentLoc)
                if value > 0.5*runningValue/iterations  
                    selectedPeaks(currentLoc) = true;
                end
            end
            break;
        end
    end
    
    % add first peak if missing
    lpFilt = lowpass(y, 2, ppcIn);
    lp = lpFilt(1:(find(selectedPeaks,1))-ppcIn/2);
    if length(lp) > 3
        [lpy, lpLoc] = findpeaks(lp);
        if ~isempty(lpy)
            [~, i] = max(lpy);
            [~, i] = min(abs(pks.l - lpLoc(i)));

            if (find(selectedPeaks,1) - i) > (ppcIn/2)
                selectedPeaks(pks.l(i)) = true;
            end

        end
    end
    
    % adjust peaks based on target number of cycles 
    for i = sum(selectedPeaks):(floor(numPeriods) - 1) % check if number of selected peaks is less than the target number, adjust w/ floor -1 for trough start/end condition
        selectedPeaksSubset.l = find(selectedPeaks); % location of all the selected peaks (indices) in the low pass filtered y, change from logical indices to linear indices
        
        values = [];
        locs = [];
        for i = 1:(length(selectedPeaksSubset.l)-1) % iterate through each selected peak (start to last, looking in between each selected peak and the next selected peak ( - 1 because there's no selected peak after the last one)
            tmp = (pks.l > selectedPeaksSubset.l(i)) & (pks.l < selectedPeaksSubset.l(i+1)); % looks for all orignally founds peaks between two adjacent selected peaks. Is a logical index of pks.l (0 if not b/n two currently selected peaks, 1 if b/n the two selected peaks currently being looked at)
            tmp = find(tmp);
            
            %left to right cost function
            locValue = 2.*abs(pks.l(tmp) - (selectedPeaksSubset.l(i) + ppcIn))/ppcIn + 1;
            locValue = 1./locValue;
            leftValue = locValue .* pks.p(tmp);
            
            % right to left cost function (location of peaks relative to 
            % target distance weighted by their prominence)
            locValue = 2.*abs(pks.l(tmp) - (selectedPeaksSubset.l(i+1) - ppcIn))/ppcIn + 1;
            locValue = 1./locValue;
            rightValue = locValue .* pks.p(tmp);
            
            % average weighting
            value = (leftValue + rightValue) / 2;
            values = [values; value];
            locs = [locs; tmp];
        end
        
        [~, maxI] = max(values);
        selectedPeaks(pks.l(locs(maxI))) = true;
        
    end
    
    selectedPeaksSubset.l = find(selectedPeaks);
    
    %% find troughs
    [~,troughs.l,~,troughs.p] = findpeaks(-y);
    selectedTroughs = zeros(size(troughs.l));
    for i = 1:(length(selectedPeaksSubset.l) - 1)
        % iterate through selected peaks. Use the cost function to find the
        % most prominent trough in between peaks
        tBetween = (troughs.l > selectedPeaksSubset.l(i)) & (troughs.l < selectedPeaksSubset.l(i+1));
        [~, index] = max(troughs.p(tBetween));
        selectedTroughs(troughs.l(find(tBetween,1)+index-1)) = true;
    end
    
    % add beginning trough if missing
    lpFilt = lowpass(-y, 2, ppcIn);
    lp = lpFilt(1:(find(selectedTroughs,1))-ppcIn/2);
    if length(lp) > 3
        [lpy, lpLoc] = findpeaks(lp);
        if ~isempty(lpy)
            [~, i] = max(lpy);
            [~, i] = min(abs(troughs.l - lpLoc(i)));

            if (find(selectedTroughs,1) - troughs.l(i)) > (ppcIn/2)
                selectedTroughs(troughs.l(i)) = true;
            end

        end
    end
    
    % add end trough if missing
    lastTrough = find(selectedTroughs, 1, 'last');
    lp = lpFilt(lastTrough:length(y));
    if length(lp) > 3
        [lpy, lpLoc] = findpeaks(lp);
        if ~isempty(lpy)
            [~, i] = max(lpy);
            [~, i] = min(abs(troughs.l - lpLoc(i) - lastTrough + 1));
            
            if (troughs.l(i) - find(selectedTroughs,1, 'last')) > (ppcIn/2 - 1)
                selectedTroughs(troughs.l(i)) = true;
            end

        end
    end
    
    selectedTroughsSubset.l = find(selectedTroughs);
    

    %% resample
    
    % critical points
    criticalPts = sort([selectedPeaksSubset.l; selectedTroughsSubset.l]);
    
    isPeak = y(criticalPts(1)) > y(criticalPts(2));
    
    
    % preallocate
    intermediateX = zeros((length(criticalPts)-1)*ppcOut/2, 1);
    currentX = 1;

    % create x values that are spaced between critical points
    for i = 1:(length(criticalPts)-1)
        if isPeak
            intermediateX(currentX:(currentX+p2t),1) = linspace(criticalPts(i),criticalPts(i+1), p2t + 1);
            currentX = currentX + p2t;
        else
            intermediateX(currentX:(currentX+t2p),1) = linspace(criticalPts(i),criticalPts(i+1), t2p + 1);
            currentX = currentX + t2p;
        end
        isPeak = ~isPeak;
    end

    % add x values for endpoints
    meanLengthp2t = 0;
    meanLengtht2p = 0;
    count = 0;
    if y(criticalPts(1)) > y(criticalPts(2))
        for i = 1:ppcOut:(length(intermediateX)-2*ppcOut)
            count = count + 1;
            meanLengthp2t = meanLengthp2t + mean(diff(intermediateX(i:(i+p2t),1)));
            meanLengtht2p = meanLengtht2p + mean(diff(intermediateX((i+p2t):(i+ppcOut),1)));
        end
        
        meanLengthp2t = meanLengthp2t/count;
        meanLengtht2p = meanLengtht2p/count;
        
        begPts = (1:meanLengtht2p:intermediateX(1,1))';
        leftValue = length(begPts);
        
        intermediateX = [begPts; intermediateX(:,1)];
        
        nanPadding = [];
        if leftValue < (peakLoc - 1)
            nanPadding = nan(peakLoc-1-leftValue,1); 
        elseif leftValue > (peakLoc-1)
            nanPadding = nan(ppcOut+peakLoc-1-leftValue,1);
        end
            
        intermediateX = [nanPadding; intermediateX(:,1)];
        
    else
        for i = 1:ppcOut:(length(intermediateX)-2*ppcOut)
            count = count + 1;
            meanLengtht2p = meanLengtht2p + mean(diff(intermediateX(i:(i+t2p),1)));
            meanLengthp2t = meanLengthp2t + mean(diff(intermediateX((i+t2p):(i+ppcOut),1)));
        end
        
        meanLengthp2t = meanLengthp2t/count;
        meanLengtht2p = meanLengtht2p/count;
        
        begPts = (1:meanLengthp2t:intermediateX(1,1))';
        leftValue = length(begPts);
        
        intermediateX = [begPts; intermediateX(:,1)];
        
        nanPadding = [];
        if leftValue < (troughLoc - 1)
            nanPadding = nan(troughLoc-1-leftValue,1); 
        elseif leftValue > (troughLoc-1)
            nanPadding = nan(ppcOut+troughLoc-1-leftValue,1);
        end
            
        intermediateX = [nanPadding; intermediateX(:,1)];
        
    end
    
    while intermediateX(end) == 0
        intermediateX = intermediateX(1:(end-1));
    end
    
    if y(criticalPts(end)) > y(criticalPts(end-1))
        intermediateX = [intermediateX(:,1); (intermediateX(end,1):meanLengthp2t:length(y))'];
    else
        intermediateX = [intermediateX(:,1); (intermediateX(end,1):meanLengtht2p:length(y))'];
    end
    
    % linearly interpolate data
    ts(:,2) = interp1(1:length(y), yInit, intermediateX);
    ts(:,1) = 0:(1/ppcOut):(length(ts(:,2))-1)/ppcOut;
    
end

function [loc,score] = chooseNextPeak(currentLoc,peakLocs,peakProm,ppc)
    i = peakLocs > currentLoc;
    peakLocs = peakLocs(i);
    peakProm = peakProm(i);
    
    locValue = 2.*abs(peakLocs - (currentLoc + ppc))/ppc + 1;
    locValue = 1./locValue;
    value = locValue .* peakProm;
    [score,loc] = max(value);
    loc = peakLocs(loc);
end