function [A,xhat,u,s,v] = phaseSpace(x,emb,tau,method)
% A = phaseSpace( x, emb, tau, (,method) )
%
% transforms the vector or matrix "x" into a phase-space representation
% governed by the choice of "method". "method" is a string that 
% can be one of the two following options:
%   "MOD" - performs Takens' method of delays (MOD) by taking time-shifted
%           versions of the columns in "x", using "emb" points spread apart
%           by "tau" samples. * "mod" is the default method
%
%   "SSA" - performs singular systems analysis (SSA) by first projecting x
%           into a higher number of dimensions (if x is a vector) or taking
%           the covariance of the columns of x (if x is a matrix), and
%           projecting down onto "emb" dimensions using SVD. If "x" is a 
%           vector, the intial embedding dimension will be 20. 
%               * optional output "xhat" is a denoised version of "x" based
%                 on the projected dimensions * 
%       
% By JMS, 10/15/2015


% check inputs
if nargin < 4
    method = 'MOD';
end

% preallocate
[n,m] = size( x );
xhat = nan;
if strcmp( method,'MOD' )
    N = n - tau*(emb-1); % maximum # of segments (k = 1,2 ... n - (e-1)*tau) 
    A = zeros( N,emb,m );
    maxE = emb;
else
    maxE = 20;
    N = n - tau*(maxE-1); % maximum # of segments (k = 1,2 ... n - (e-1)*tau)
    A = zeros( N,maxE,m );
end 

% loop through series "x" and pull out e-dimensional vectors
for e = 1:maxE
    startInd = tau*(e-1) + 1;
    stopInd = tau*(e-1) + N;
    A(:,e,:) = x(startInd:stopInd,:);
    %A(:,e,:) = [x(stopInd+1:end),x(1:startInd-1,:);x(startInd:stopInd,:)]; % indexes by tau and allows for circular modes 
end

% check if "SSA", and perform PCA to reconstruct "emb" principal dimensions
if strcmp(method,'SSA') 
    A = A / (size( A,1 )^.5); % normalize
    A = reshape( A,size( A,1 ),m*maxE );
    [u,s,v] = svd( A,'econ' );
    
    % project the data onto lower dims using the PC's in "v" and the original A matrix
    A = u(:,1:emb) * s(1:emb,1:emb); % u*s...the projected data onto time
    
    % form denoised reconstruction xhat using reduced dims and solving usv'
    if nargout > 1
        xhat = mean( A * v(1:emb,1:emb)',2 ); % u*s*v'
        xhat = xhat * size( xhat,1 )^0.5;
    end
    
    % print the corresponding variance explained
    varExp = cumsum(diag(s)/sum(diag(s)));
    fprintf('projected onto %0.3f percent of variance\n',varExp(emb)*100);
end

end