function D = phaseSpaceDist( A,varargin )
% D = phaseSpaceDist( A, (A2) )
% 
% compute the euclidean distances between every vector pair in the phase
% space matrix A. You can optionally provide a separate phase space matrix
% A2, in which case D will be computed via the distances between vector
% pairs of A and A2 (known as the cross-distance matric, to compute the
% cross-recurrence)
%
% D is then normalized to contain values between [0,1] via D / max( D )
  
% check for A2
if nargin < 2 || isempty(varargin{1})
    A2 = A;
else 
    A2 = varargin{1};
end

% compute pairwise distances
D = compute_pairwise_dist( A,A2 );

% normalize the distance matrix
D = D / max( D(:) );

end