function [fd,gf,L] = hfractal(X,varargin)
% ----------------- [fd,gf,L] = hfractal(X,varargin) ------------------
%
%   Calculates the fractal dimension of input matrix/vector "X" using
%   Higuchi's algorithm (i.e. via time-lagged vectors from each
%   column of X). This is a vectorized version to increase 
%   speed, especially pertinent when K or X are large.
%
% Inputs:
%   X   : column-oriented vector or matrix of data
%   K    : (optional) the maximum iterative value for time-delay
%            reconstruction...default = min(100,size(X,1)/2);
%   pl   : (optional) plotting as 1 or 0...default = 0;
%
% Outputs:
%   fd  : the fractal dimension of each column in X
%   gf  : goodness of fit of least-squares line to lengths in L
%   L   : the distances calculated by the higuchi algorithm. Rows are
%          lengths for k = 1:Kmax, columns are for each time seires in X
%
% By JMS, 4/6/2016 
%-------------------------------------------------------------------------

% check inputs
if max(size(X)) == 1 && isrow(X)
    X = X';
end
if nargin > 1 && ~isempty(varargin{1})
    K = varargin{1};
else
    K = floor(min(200,size(X,1)/2)); 
end
if nargin > 2 && ~isempty(varargin{2})
    pl = varargin{2};
else
    pl = 0; 
end

% get sizes and preallocate
N = size(X,1);
M = size(X,2);
L = zeros(K,M);
fd = zeros(1,M);

xpoints = log(1:K)'; % for getting slope of lengths
gf = zeros(1,M); % for goodness of fit 

% begin iterative process
for i = 1:M
    for k = 1:K
        m = 1:k;
        lim = min(fix((N-m) / k)); % get minimum limit to avoid indexing out of bounds
        ind = repmat([0:k:lim*k]',1,k);
        ind = bsxfun(@plus,ind,m);
         
        % make matrix of this column of X, expanded "k" times
        lag = repmat(X(:,i),1,k); 
        lag = lag(ind); % matrix of time-lagged vectors 
        
        % subtract xlag_i-1 from xlag_i and sum to obtain length
        normalize = (N - 1) / (lim * k);
        dis = sum( abs(lag(2:end,:) - lag(1:end-1,:)) ); % sum( | x(m + ik) - x(m + (i-1)k) | )
        dis = dis*normalize / k; % 1/k * dis * (N-1) / (N-m/k)*k
       
        % store the mean of subseries m = 1...k into index k of L
        L(k,i) = mean(dis);
        
        clear dis xlag ind lim normalize
    end
    
    % now calculate the least-squares slope of the line spanned by L for
    % time series i
    ypoints = log(L(:,i));
    f = polyfit(xpoints,ypoints,1); % y = mx + b
    fd(i) = -f(1); % slope
    yhat = f(1)*xpoints + f(2);
     
    % optional output...goodness of fit
    if nargout >=2
        gf(i) = goodnessOfFit(yhat,ypoints,'nrmse');
    end
    
    % optional plotting
    if pl == 1
        figure; hold on
        scatter(xpoints,ypoints);
        plot(xpoints,yhat,'r','linewidth',2);
        title(['fractal dim: ',num2str(-f(1))]);
    end
    
    clear ypoints f yhat

end % 1:M (columns of X)

end
   