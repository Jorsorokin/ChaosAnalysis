function [mH,H,P,C] = mmPerm(X,tau,m,varargin)
% -------------- [mH,H,P,C] = mmPerm(y,tau,m,varargin) -------------------
%
%   Calculates the multi-variate (and optionally multiscale) permutation
%   entropy from a matrix of column-oriented vectors X using the approach
%   outlined in Morabito et al. 2012. 
%
%   Multi-variate permutation entropy assumes that multiple channels in X
%   are not independent and together provide information regarding the
%   complexity of each individual signal X.
%
%                   >>> INPUTS >>>
%
% X     : the data matrix (column-oriented)
% tau  : the time lag
% m    : the embedding dimension
% S     : OPTIONAL...the scales to use. If provided, multi-scale
%           permutation entropy will be calculated for each column in X.
% normalize : OPTIONAL...will normalize entropy by log2(m!) (default = 0)
%
%                   <<< OUTPUTS <<<
%
% mH : the multi-variate permutation entropy
% H    : permutation entropy for each column in X over multiple scales  
% P     : the probability matrix of each motif 1...m! in each column
%          divided by the number of columns (thus sum(sum(P)) = 1)
% C     : contingency as the MSE between mH and mean(H)
%
% By JMS, 4/8/2016
%------------------------------------------------------------------------------------

% check inputs
if nargin > 3 && ~isempty(varargin{1})
    S = varargin{1}; 
else S = 1; end % if S = 1, same result as non multi-scale permutation entropy
if nargin >4 && ~isempty(varargin{2})
    normalize = varargin{2};
else normalize = 0; end

% presets/scalars
N = size(X,2); % number of columns
H = zeros(S,N);
mH = zeros(S,1);
P = zeros(factorial(m),N,S); % m! x N x S matrix

for s = 1:S
    % resample by scale s
    rX = course_grain(X,s);
    
    % permutation entropy
    [H(s,:),p] = perm_entropy(rX,tau,m);
    
    % divide probabilities by num columns C
    P(:,:,s) = p / N; % sum(sum(P)) == 1
    
    % get multi-variate permutation entropy
    mP = sum(P(:,:,s),2); % sum across channels
    mH(s,1) = -sum(mP .* log2(mP+eps));
    
    clear p mP
end

% normalization
if normalize == 1
    H = H / log2(factorial(m));
    mH = mH / log2(factorial(m));
end

% contingency
if nargout == 4
    mH2 = mean(H,2);
    C = mean( (mH - mH2).^2 ); % mse between mH and mean(H)
end

end

%% Functions

% Permutation entropy
%===========================
function [H,p,count] = perm_entropy(D,tau,m)
    % calculate permutation entropy for all columns in X
    %
    % Inputs:
    %   D: data matrix
    %   tau: time delay
    %   m: embedding dimension
    %
    % Outputs:
    %   H: permutation entropy
    %   p: probabilities of each motif (row) for each column in D
    %           h / (N - m + 1)
    %   count: raw counts of each motif for each column in D

    % get limit for iterating, create empty list for storing permutation
    % entropy of each column
    N = size(D,1);
    C = size(D,2);
    lim = N - (m-1)*tau;
    plist = perms(1:m)'; % rows are motifs, columns are 1...m!
    L = size(plist,2);
    count = zeros(L,C); % will store the frequencies of occurence of each motif for each column
    
    % loop through and get the frequencies of each motif
    for i = 1:lim
        [~,motif]=sort(D(i:tau:[i + tau*(m-1)],:));
        for v = 1:C
            check = sum( abs(bsxfun(@minus,motif(:,v),plist)) );
            count(:,v) = count(:,v) + (check' == 0);
        end
    end
    
    % get probability of the counts
    p = bsxfun(@rdivide,count,sum(count));
    
    % calculate permutation entropy H
    H = -sum(p .* log2(p+eps)); % add "eps" to avoid log2(0)...doesn't change actual calculations much
end
%===========================


% course multi-scale series
%===========================
function rD = course_grain(D,S)
    % resample vectors in D by 1/S

    N = length(D);
    J = fix(N/S);
    rD = zeros(floor(N/J),size(D,2));

    for i=1:J
        rD(i,:) = mean(D((i-1)*S+1:i*S,:),1);
    end

end
%===========================