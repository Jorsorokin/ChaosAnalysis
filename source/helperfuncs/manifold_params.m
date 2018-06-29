function [emb,tau,gammaBar,S] = manifold_params( x )
% [emb,tau,gammaBar,S] = manifold_params( x )
%
% uses the average displacement (AD) and gamma test to estimate the best 
% tau and embedding dimensions from the time series in "x"
%
%               >>> INPUTS >>>
% x: 
%   column-oriented time series vector
% 
%               <<< OUTPUTS <<<
% emb:
%   the estimated embedding dimension
% tau:
%   the estimated time delay 
% gammaBar:
%   the y-intercept from the gamma test...where the differences between
%   gammaBar decreases to 10% is the indicated embedding dimension
% S:
%   the results from the average displacement test...for each column, where
%   the differences between S decrease to 40% is the best tau
%
% By JMS, 12/15/2016

% loop over embedding and tau to estimate both
embedding = 3:30;
nEmb = numel(embedding);
gammaBar = zeros(nEmb,1); % 16 - 3 embedding dims
s = zeros(nEmb,1);
tau = zeros(nEmb,1);
nTau = 20;
S = zeros(nTau,nEmb);

for e = 1:numel(embedding)
    m = embedding(e);
    
    %% Part 1: estimate tau 
    for t = 1:nTau
        A = phaseSpace(x,m,t);
        dA = bsxfun(@minus,A(:,2:end),A(:,1));
        dA = (sum(dA.^2,2)).^0.5;
        S(t,e) = mean(dA);
    end
    
    % now find optimal tau given m
    slope = diff(S(:,e))./max(diff(S(:,e)));
    tau(e) = find(slope <= .4,1); % where slope < .4 of max slope

    %% Part 2: perform the gamma test so that we can estimate best m and tau
    eta = phaseSpace(x,m,tau(e)); % the input, given the best estimated tau
    M = size(eta,1);
    y = eta(:,1)+1; % the output
    maxP = min(M,20);
    etaDist = zeros(maxP,1);
    yDist = etaDist;
    
    % get distances and 20-nearest neighbors
    for i = 1:M
        d = bsxfun(@minus,eta,eta(i,:));
        d = (sum(d.^2,2)).^.5; % euclidean distances
        d(i) = 1000000; % arbitrarily large
        [~,loc] = sort(d); % sort the distances
        etaDist = etaDist + d(loc(1:maxP)).^2; % add the squared distances
        yDist = yDist + (y(loc(1:maxP)) - y(i)).^2; % ditto for output
    end
    etaDist = etaDist / i; % divide by i to take the mean
    yDist = yDist / (i*2); % ditto
    
    % now take the running averages across p
    delta = cumsum(etaDist) ./ (1:maxP)';
    gamma = cumsum(yDist) ./ (1:maxP)';
    
    % now fit a line to (p,gamma) to get the slope and yint
    B = delta \ [ones(size(gamma)),gamma];
    s(e) = B(2);
    gammaBar(e) = B(1);
end

%% Part 3: get m and tau 
dGamma = diff(gammaBar)/min(diff(gammaBar));
emb = find(dGamma <= .1,1);
tau = tau(emb);

end


