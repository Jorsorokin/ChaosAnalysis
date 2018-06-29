function leField = FTLE( A,varargin )
    % Compute the finite time lyapunov exponent (FTLE) field
    % 
    % The maximum lyapunov exponent for a given point / trajectory
    % describes the exponential growth of nearby trajectories in a flow field.
    % Lyapunov exponents < 1 represent shrinking growth (attracting
    % trajectories) and > 1 describes unstable positive growth (diverging
    % trajectories). By computing the lyapunov exponent for each point in time,
    % one can create a lyapunov field for a time-evolving trajectory that
    % relates different regions of underlying manifold from which the
    % particular state-space matrix A is realized with different types of
    % trajectory behaviors.
    %
    % Usage:
    %   leField = FTLE( A,(radii,nT,minT) )
    %
    % Inputs:
    %   A - phase space matrix
    %
    %   (radii) - the radius (in percentage) of the cloud surrounding each state
    %             vector in A, which determines the nearest points for
    %             calculating the transition matricies (default = 5% of range of
    %             each dimension)
    %
    %   (nT) - the number of time steps to evolve the point cloud (default = 5)
    %
    %   (minT) - minimum time steps away from reference point to consider another
    %          point in the cloud. Helps avoid autocorrelation (default = first min of autocorrelation)
    %
    %   (dt) - the number of time steps to iterate by. Default = 1
    %
    % Outputs:
    %   leField - layapunov exponent field computed over the 
    %             phase space trajectory A. Large values of the leField
    %             describe repelling regions in the phase-space that separate
    %             different behaviors of the trajectory
    % 
    % By JMS, 1/3/2017

    % get some parameters of the phase space topology
    [nPoints,nDim] = size( A );
    r = range( A );

    % check the inputs
    if nargin > 1 && ~isempty(varargin{1})
        radii = varargin{1};
    else 
        radii = 0.05;
    end

    if nargin > 2 && ~isempty(varargin{2})
        nT = varargin{2};
    else
        nT = 5;
    end

    if nargin > 3 && ~isempty( varargin{3} )  
        minT = varargin{3};
    else
        % find the first minimum of the autocorrelaton to find the number of time
        % steps to skip surrounding each point when finding nearest neighbors
        xA = xcorr(A(:,1)); % use the first column in case of SSA method, in which 
                            % most of the variance will be contained in the first
        xA = xA(nPoints:end);
        minT = find( max(xA) ./ xA <= .3,1 );
    end

    if nargin > 4 && ~isempty( varargin{4} )
        dt = varargin{4};
    else
        dt = 1;
    end

    % expand the radii into a matrix to develop the point cloud
    radii = r * radii;
    pointMat = zeros( 2*nDim,nDim );
    for j = 1:nDim
        pointMat(j*2-1:j*2,j) = [radii(j);-radii(j)];
    end

    % loop over the state vectors in A to form the layapunov exponent field
    leField = nan( 1,nPoints-nT ); 
    
    % indicies to loop
    start = 1;
    stop = nPoints - nT;
    
    for p = start:dt:stop

        % find the nearest points to each of the vertices
        vertices = bsxfun( @plus,A(p,:),pointMat ); % encapsulating cloud
        removeInds = [max(1,p-minT),min(nPoints-nT,p+minT)];
        
        % compute distances
        dist = compute_dist( A,vertices );
        dist([removeInds(1):removeInds(end),nPoints-nT+1:nPoints],:) = 1e6; % arbitrarily large

        for j = 1:nDim*2
            [~,loc(j)] = min(dist(:,j)); % avoids redundant neighbors
            dist(loc,j+1:end) = 1e6;
        end

        % now get the distances of the nearest neigbors to each other
        % at time t0 and at time t0 + nT 
        neighbors_t0 = A(loc,:); % nearest point to each vertex that is "minT" points away of the reference point
        neighbors_t = A(loc+j,:);

        % compute the jacobian of the deformation tensor
        d_t0 = diag( neighbors_t0(1:2:end,:) - neighbors_t0(2:2:end,:) );
        d_t = neighbors_t(1:2:end,:) - neighbors_t(2:2:end,:);
        dF = (d_t ./ d_t0)';

        % now find the spectral norm of the deformation matrix, which
        % describes the maximum growth rate of the local vector field.
        s = svd( dF,'econ' ); % equivalent to: sqrt( eig( dF'dF ) ) - eigendecomp of the Cauchy-Green Strain Tensor
        leField(p) = s(1);
    end

    % finally, compute the lyapunov exponents from the eigenvalues
    leField = 1/nT * log( leField )'; 
   
    
    %% Helper function
    function D = compute_dist( X,Y )
        D = real( sqrt( bsxfun( @plus, sum( X.^2,2 ),...
            bsxfun( @minus, sum( Y.^2,2 )', 2*(X * Y') ) ) ) );
    end
    
end

