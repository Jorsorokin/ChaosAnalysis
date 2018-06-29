function [CR,D] = computeCrossRecurrence(A1,A2,radius)
% [CR,D] = computeCrossRecurrence( A1,A2 );
% 
% find the distances between vectors of the phase space matrices A1 and A2
% to compute D, then take the cross recurrence by thresholding D. Can
% optionally specify "radius" between [0,1] to threshold D
%      
% D is then normalized to contain values between [0,1] via D / max(max(D))
  
if nargin < 3
    radius = [];
end

% preallocate
N = size(A1,1);
M = size(A2,1);
D = zeros(N,M); 

% loop through A and get distances with other vectors in A
for i = 1:M
    distances = bsxfun(@minus,A1,A2(i,:)); % distance between vector A2(i) and all of A1
    D(:,i) = sqrt(sum( distances.^2,2 )); % euclidean distance
end

% normalize to max to make all elements in [0, 1] form
D = D / max(max(D));

% now compute the cross recurrence
CR = computeRecurrence(D,radius);

end