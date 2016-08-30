function P = plot_circle(x,y,r,color)
% Plot a circular arc as a pie wedge.
% a is start of arc in radians, 
% b is end of arc in radians, 
% (h,k) is the center of the circle.
% r is the radius.
% Try this:   plot_arc(pi/4,3*pi/4,9,-4,3)
% Author:  Matt Fig
t = 0:0.05:2*pi;
xp = r*cos(t) + x;
yp = r*sin(t) + y;
%x = [x h x(1)];
%y = [y k y(1)];
P = fill(xp,yp,color); hold on;
if ~nargout
    clear P
end

end