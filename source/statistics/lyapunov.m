function [le,gf] = lyapunov(A,D,Fs,varargin)
% --------------- [le,gf] = lyapunov(A,D,Fs,varargin) ---------------
%
%   Calculates the largest lyapunov exponent "le" from the state-space
%   matrix "A" and the distance matrix "D" (i.e. the matrix of distances
%   between state-space vectors in A).
%
%   Uses Rosenstein's algorithm for determining the maximum lyapunov exp. 
%
% Inputs:
%   A   : state-space matrix (column-oriented)
%   D   : distance matrix 
%   Fs  : the sampling rate
%   md  : OPTIONAL...minimum distance or time series vector. 
%         If vector, will calculate the minimum distance as the 
%         reciprocal of the mean frequency of the fourier transform
%   pl  : plotting...1 or 0 (default = 0); 
%
% Outputs:
%   le   : the largest lyapunov exponent
%   gf   : goodness-of-fit of equation y = mx + b 
%
% By JMS, 4/5/2016
%-----------------------------------------------------------------------------

% optionals
if nargin > 3 && ~isempty(varargin{1})
    md = varargin{1};
else md = .01; end
if nargin > 4 && ~isempty(varargin{2})
    pl = varargin{2};
else pl = 0; end

% indexing variables
M = size(A,1);
maxiter = min(400,M/2); % 400 should be enough points for least squares slope
distance = zeros(maxiter,1);


% get min distance via 1/mean frequency of "x" if "d" is a vector
if numel(md) > 1
    if isrow(md)
        md = md';
    end
    npoints = numel(md);
    fourier = fft(md); % take the fourier trasnform to extract power for different frequenceis
    fourier = fourier(1:npoints/2+1,:); % only take positive frequency components
    freqs = 0:Fs/npoints:Fs/2; % frequency vector based on the sampling rate
    mf = freqs * log10(abs(fourier)) ./ sum(log10(abs(fourier))); % the mean frequency
    md = 1/mf;
    clear fourier freqs
end


% normalize distance matrix
D = D / max( D(:) );

% make main diagonal from distance matrix large (i.e. 1000)
I = eye(M) * 1000;
D = D + I;
clear I

% get minimum distance/position of each column D
% after excluding distances less than "md" 
D(D<md) = 1000; % make arbitrarily large
[~,pos] = min(D);


% loop through A and compare to min distances of D
for k = 1:maxiter
    maxind = M - k;
    p = pos(1:maxind); % position vectors of min(D)
    invalid = p>maxind; % any position vectors > max index
    valid = p(p<=maxind); % any positions <= max index
    a1 = A([1:maxind]+k,:); % j = 1,2...maxind
    a2 = A(valid+k,:); % only position vectors satisfying maxind
    a1(invalid,:) = []; % remove bad indexes from a1
    
    % distances between a1,a2 vectors 
    dis = sqrt(sum( (a1 - a2).^2, 2)); 
    dis = dis(dis>0);
    
    if ~isempty(dis)
        distance(k) = mean(log(dis));
    else
        distance(k) = 0;
    end
    
    clear a1 a2 p valid invalid
end


% calculate the maximum lyapunov exponent via least-squares fit of "distance"
xpoints = 15:maxiter-15;
time = xpoints/Fs;
params = polyfit(time',distance(xpoints),1); % y = mx + b
le = params(1);
y = params(1)*time' + params(2);

% goodness of fit
if nargout == 2
    gf = goodnessOfFit(distance(xpoints),y,'nrmse');
end

% optional plotting of distance and fit
if pl == 1
    figure; plot((1:numel(distance))/Fs,distance);
    hold on; plot(time,params(1)*time + params(2),'r','linewidth',2);
    title(['Lyapunov exponent: ',num2str(le)]);
end


end

