function Q = RQA(RP,minL)
% =========================================================================
% Q = RQA( RP (,minL) )
%
% Quantification of recurrence plots for nonlinear dynamical analysis
%
% INPUTS:
%   RP : the symmetric recurrence plot
%
%   minL : OPTIONAL...the minimum length of diagonal/vertical line
%          for quantification. Default = 2
%
% OUTPUTS:
%   Q: a structure with the following fields:
%       N:
%           the max diagonal size of the RP
%       dCount:
%           frequency dist. of lengths of diagonals.
%       vCount:
%           frequency dist. of lengths of verticals
%       Pd:
%           probability dist. of "diagCount"
%       Pv:
%           probability dist. of "vertCount"
%       RR:
%           recurrence rate - percentage of "1" in the RP
%       det:
%           determinism - percentage of points in diagonals
%       lam:
%           laminarity - ratio points in vertical lines vs. single points
%       TT:
%           trapping time - average length of vertical lines
%       L:
%           length - average length of diagonal lines
%       div:
%           divergence - 1 / max(diagonal line length)
%       entr:
%           entropy - shannon entropy of "Pd"
%
% By JMS, 10/31/16
% =========================================================================

%% Params/presets
if nargin < 2
    minL = 2; % minimum length of diagonal/vertical for calculations
end
N = size(RP,1); % the length of dim. of the RP plot
diagCount = zeros(N,1); % for storing freq distribution of diagonal lengths
diagIDX = (1:N)'; % actual diagonal lengths possible
vertCount = zeros(N,1); % for storing freq distribution of vertical lengths
vertIDX = diagIDX; % actual vertical lengths possible

%% loop over diags/columns of RP
for k = 1:N

    % extract the diagonal and vertical values
    diagonal = [0; diag( RP,k ); 0]; % extract the diagonal elements of the kth diagonal
    vertical = [0; RP(1:k,k); 0]; % extract the vertical, up to the major diagonal

    % get the counts of the line segments in each
    if max( diagonal ) > 0
        diagCount = get_lengths( diagonal,diagCount );
    end
    if max( vertical ) > 0
        vertCount = get_lengths( vertical,vertCount );
    end
end

% multiply the counts by 2, since we only quantified half of the recurrence
% plot (since it is symetric)
diagCount = diagCount * 2;
vertCount = vertCount * 2;

%% PERFORM RQA ANALYSIS

% set presets for the variables, which will change if sr > 0
det = 0;
lam = 0;
L = 0;
tt = 0;
Pd = 0;
Pv = 0;
entr = 0;
div = 0;

% ==================================
% sum of recurrence (sr) & recurrence rate (rr)
% ==================================
sr = sum( sum( RP ) );
rr = sr / N^2;

if sr ~= 0
    % ==================================
    % determinism (det)
    % ==================================
    totalDiag = diagCount .* diagIDX;
    det = sum( totalDiag(minL:end) ) / sr;

    % ==================================
    % laminarity
    % ==================================
    totalVert = vertCount .* vertIDX;
    lam = sum( totalVert(minL:end) ) / sum( totalVert(1) );

    % ==================================
    % probabilities (Pv & Pd)
    % ==================================
    Pd = diagCount / sum( diagCount );
    Pv = vertCount / sum( vertCount );

    % ==================================
    % entropy (entr)
    % ==================================
    entr = -nansum( Pd(minL:end) .* log( Pd(minL:end) ) ); % shannon entropy

    % ==================================
    % Averaged diagonal length & trapping time (avg. vert. len)
    % (L & TT)
    % ==================================
    L = sum( totalDiag(minL:end) ) / sum( diagCount(minL:end) );
    tt = sum( totalVert(minL:end) ) / sum( vertCount(minL:end) );

    % ==================================
    % divergence (div) as lyapunov approximation
    % ==================================
    Lmax = max( diagIDX(diagCount > 0) );
    div = 1/Lmax;
end

% Finally store everything into the structure Q
Q = struct('N',N,'minL',minL,'dCount',diagCount,'vCount',vertCount,...
           'Pd',Pd,'Pv',Pv,'RR',rr,'det',det,'lam',lam,'TT',tt,...
           'L',L,'div',div,'entr',entr);

end

%% Helper functions

% get_lengths()
function count = get_lengths( line,count )
    % get the lengths of lines-segments, and their frequency of occurrence
    % in a vertical/diagonal segment of the recurrence matrix

    % differentiate to find lengths of diagonal lines
    lineDiff = [0; diff( line ); 0];
    idx = 1:numel( lineDiff );

    % check the elements of lineDiff...find lengths
    if max( lineDiff ) == 0 % if 1 long diagonal, set L = to the sum of the points
        lineLength = sum( line );
    else
        on = idx(lineDiff == 1);
        off = idx(lineDiff == -1);
        lineLength = off - on; % find the lengths of the diagonal elements
    end

    % count the number of times a segment of length "d" occurred.
    % store into "count" to use for quantification of DET, etc...
    uL = unique( lineLength );
    for d = 1:numel( unique( lineLength ) )
        count(uL(d)) = count(uL(d)) + sum( lineLength == uL(d) );
    end
end
