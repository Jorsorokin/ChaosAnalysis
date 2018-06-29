function [sH,lH] = plotLCS( X,leField )
% [sH,lH] = plotLCS( X,leField )
%
% Plot the lyapunov field and lagrangian coherent structures calculated by
% "LCS.m" onto a gridded, interpolation of the bounds of X. If X has > 2
% dimensions, plots onto the first two dimensions

n = size( X,1 );
m = numel( leField );
if m < n
    X = X(1:m,:);
end
g = fspecial( 'gaussian',21,3 );
bounds = [min( X );max( X )];
[xq,yq] = meshgrid( linspace( bounds(1,1),bounds(2,1),250 ),...
                    linspace( bounds(1,2),bounds(2,2),250 ) );
Vq = griddata( X(:,1),X(:,2),leField,xq,yq );
Vq = imfilter( Vq,g );

sH = surf( xq,yq,Vq,'linestyle','none' );
hold on
lH = plot3( X(:,1),X(:,2),leField+0.1,'k.-' );

end