function [state,tvec] = recurrence_movie(x,Fs,window,overlap,tau,emb,mov)
% [state,time] = recurrence_movie( x,Fs,window,overlap,tau,emb,(,mov) );
%
% calculates recurrence parameters using RQA on small blocks of the data in
% "x" using "overlap" amount of points for the sliding window. If "mov" is
% left blank, will automatically display the results of the sliding window,
% else if "mov" is 0, will not plot.
%
% returns an N x 6 matrix in "state", with rows equal to the number of
% split data segments, and columns equal to:
%   | det | L | lam | TT | entr | div |
% 
% also returns the time vector "time", which contains the middle time point
% of the sliding window from which the variables in "state" were derived

% check mov input
if nargin < 7 || isempty(mov)
    mov = 1;
end

window = round(window*Fs);
overlap = round(overlap*Fs);
[bounds,xsplit] = splitdata_overlap(x,window,overlap);

% initialize/preallocate
[n,~] = size(x);
nvec = size(xsplit,2);
time = linspace(0,n/Fs,n);
tvec = zeros(nvec,1);
div = tvec;
det = tvec;
TT = tvec;
L = tvec; 
entr = tvec;
lam = tvec;
Alim = [-5 5];

% loop
for j = 1:nvec
    tvec(j) = mean(bounds(j,:))/Fs;
    
    % get the phase-space representation
    A = phaseSpace( xsplit(:,j),emb,tau );
    
    % perform SVD and project A to lowest 3 dimensions in order to 
    % project the trajectory into a space with the highest variance
    [u,s,~] = svd( A );

    % project down (U'S) to get at least 50% variance explained
    sigma = diag( s );
    sigma = cumsum(sigma)/sum(sigma); % gets a vector of variances explained
    ndim = find(sigma > 0.5, 1);
    if ndim < 3
        ndim = 3;
    end
    A = u(:,1:ndim) * s(1:ndim,1:ndim);
    
    % get the recurrence matrix and quantification
    D = phaseSpaceDist( A' );
    R = computeRecurrence( D );
    Q = RQA( R,3 );
    div(j) = Q.div;
    det(j) = Q.det;
    TT(j) = Q.TT;
    L(j) = Q.L;
    entr(j) = Q.entr;
    lam(j) = Q.lam;
    
    % decide if plotting
    if mov == 1
        % plot the data and the sliding window
        subplot(2,1,1); 
        timeWin = time(bounds(j,1):bounds(j,2));
        plot(time,x,'k');hold on
        plot([timeWin(1) timeWin(1)],[-8 4],'r',[timeWin(end) timeWin(end)],[-8 4],'r');
        set(gca,'ylim',[-8 4],'xlim',[time(1) time(end)]);
        hold off

        % plot the phase-space trajectory
        subplot(2,3,4);
        plot3(A(:,1),A(:,2),A(:,3),'k');
        set(gca,'xlim',Alim,'ylim',Alim,'zlim',Alim);
        title(sprintf('D_{proj}: %02d',ndim));

        % plot the recurrence matrix 
        subplot(2,3,5);
        imagesc(1-R); colormap('gray');
        set(gca,'yticklabel',[],'xticklabel',[]);

        % plot the divergence "div"
        subplot(4,3,9);
        plot(tvec(1:j),div(1:j),'ko-'); hold on
        plot(tvec(1:j),smooth(div(1:j)),'r','linewidth',2); hold off

        % plot the entropy "entr"
        subplot(4,3,12);
        plot(tvec(1:j),entr(1:j),'ko-'); hold on
        plot(tvec(1:j),smooth(entr(1:j)),'r','linewidth',2); hold off;

        drawnow
    end
end

if mov == 1
    suptitle(sprintf('window: %02d, overlap: %02d, tau: %02d, emb: %02d',window,overlap,tau,emb));
end

% store into the "state" matrix
state = horzcat(det,L,lam,TT,entr,div);

end