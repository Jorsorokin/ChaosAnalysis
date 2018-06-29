function G = RGA(RP)
% =========================================================================
% G = RGA( RP )
%
% Graph-network analysis of recurrence plots (see Donner et al. 2010)
% 
% extracts quantitative variables from the recurrence plot using
% graph-theoretical approaches. Each index i in the recurrence
% plot is considered a "vertex", and other indices j, where i ~= j, are 
% connected "nodes" if that point [i,j] in the recurrence plot is 1. Using
% this framework, the following ouputs are calculated, which complement
% traditional recurrence quantification analysis (RQA) by providing extra
% measures on the local and global topological features of the phase-space
% trajectory used to construct the RP
%
% INPUTS:
%   RP : the symmetric recurrence plot
%
%   minL : OPTIONAL...the minimum length of diagonal/vertical line
%          for quantification. Default = 2
%
% OUTPUTS:
%   G: a structure with the following fields:
%       k:
%           SUM{RP_v} for column v (sum of all recurrence points to j)
%       rr:
%           the local recurrence rate for each vertex (also known as the 
%           "degree centrality" of the vertex v), defined as 1/N-1 *
%           SUM{RP_i} = 1/N-1 * kv, where N = the number of rows
%       ent:
%           the shannon entropy of the vertex degrees "k(i)" defined as:
%               -SUM{P(k) * log{P(k)}}
%       c:
%           the local clustering coefficient, defined as the number of
%           nodes connected to vertex "v" that are also connected to
%           themselves (i.e. the number of closed triangles including
%           vertex "v"). Cv is between [0 1], defined by the equation: 
%               c(v) = (2*nt(v)) / (kv(kv-1)), where Nv = # of closed triangles
%       nt:
%           the number of connected triangles for each vertex, used to
%           calculate the clustering coefficient "c"
%       l:
%           characteristic path length
%       e:
%           characteristic efficiency
%
% By JMS, 10/31/16
% =========================================================================

% preallocate
[N,M] = size( RP );

% calculate k and rr, since they don't need a loop
k = sum( RP ); % subtract 1 to eliminate i=j recurrence

% calculate nt(v) connected triangles
nt = sum(RP^3 .* eye(N));

% now compute cluster coefficient
c = nt ./ (k .* (k-1));
c(isinf(c)) = 0;

% calculate rr
k = k-1; % eliminate i==j
rr = k / (N-1);

% calculate entropy
Pk = k/sum(k);
ent = -nansum(Pk.*log(Pk));

% shortest path distances between each node
D = graphallshortestpaths(sparse(RP));

% characteristic path length
l = nansum(D(:))/(N*(N-1));

% global efficiency
e = (nansum(nansum(1./(D+eye(N)))) - N)/(N*(N-1));

% store
G.k = k;
G.rr = rr;
G.nt = nt;
G.c = c;
G.ent = ent;
G.l = l;
G.e = e;

end
        