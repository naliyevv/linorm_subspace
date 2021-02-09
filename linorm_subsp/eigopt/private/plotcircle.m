function [] = plotcircle(rad,center,plotfig)
% call : plotcircle(rad,center,plotfig)
%   task
%                 Draws circle with the given radius.
%
%   input
%       rad     : radius of the circle which will be drawn
%       center  : center of the circle  
%       plotfig : figure handle, if no figure with this handle
%                 exists, a new figure will be opened, otherwise
%                 the circle will be drawn on the figure with
%                 handle plotfig.
%   output
%                 no output

num_of_points = 1000; 
rads = rad * rad;
xv = [];
yv = [];
inc = 4*rad / num_of_points;

rx = real(center);
ry = imag(center);



x = -rad;
while x < +rad
    xv = [xv x+rx];
    y = sqrt(rads - x*x);
    yv = [yv y+ry];
    x = x + inc;
end

x = +rad;
while x > -rad
    xv = [xv x+rx];
    y =  -sqrt(rads - x*x);
    yv = [yv y+ry];
    
    x = x - inc;
end

xv = [xv xv(1)];
yv = [yv yv(1)];


figure(plotfig);
plot(xv,yv,'b-');