function R = computeRecurrenceNeighbor( D,nNeighbor )
% R = computeRecurrenceNeigbor( D, (nNeigbor) );
% 
% computes the recurrence matrix R from the distance matrix D by finding
% a fixed number of nearest neighbors around each point in D (locally
% constrained recurrence rate). If "nNeighbor" is left blank, will
% automatically find 20 nearest neighbors

% check inputs
if nargin < 2 || isempty(nNeighbor)
    nNeighbor = 20;
end

% presets
N = size( D,1 ); 
indices = reshape( 1:N*N,N,N ); % for extracting diagonal elements of D
R = zeros( N,N );

% clean up D
for j = 1:3 % consecutive time points
    D(diag( indices,j )) = 1e6;
    D(diag( indices,-j )) = 1e6;
end

% loop over points in D and find nearest "nNeighbor" other points
for p = 1:N
    [~,loc] = sort( D(:,p) ); % sort ascending order
    R(loc(1:nNeighbor),p) = 1;
end

% convert to logical to save space
R = logical( R );

end