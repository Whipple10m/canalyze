function handle = arrow(start,stop,scale)

%  ARROW(start,stop,scale)
%                          
%  Draw a line with an  arrow at the end of a line
%  start is the x,y point where the line starts
%  stop is the x,y point where the line stops
%  Scale is an optional argument that will scale the size of the arrows
%  It is assumed that the axis limits are already set

%       9/5/96    Steve Wineberg
%       1. Modified from arrow2.m to scale arrow heads appropriately in
%          x and y if xd does not equal yd.
%       2. Modified to make arrow slope consistent with xd and yd scaling
%       3. Also modified to avoid possible divide by xd = 0.


% get axis ranges set before call to this routine.
  xl = get(gca,'xlim');
  yl = get(gca,'ylim');
  xd = xl(2)-xl(1); 
  yd = yl(2)-yl(1); 

% this sets the scale for the arrow size, thus enabling the arrow 
% to appear in correct proportion to the current axis
xdif = stop(1) - start(1);
ydif = stop(2) - start(2);
if nargin==2
  scale = .006;    % Or whatever value gives desired arrowhead size
end

axis(axis)

if( (xdif == 0) & (ydif < 0) )       % In case xdif is zero.
  theta=-pi/2;
elseif( (xdif == 0) & (ydif > 0) )
  theta=pi/2;
else
  theta = atan(ydif/xdif);  % actual slope of arrow on the graph 
                            % depends on xd, yd, as well as xdif, ydif
end

if(xdif>=0)
  scale = -scale;
end

% xd and yd are x and y axis expansion factors, respectively.
% multiply x-component of arrowhead by xd and y-component by yd so 
% arrowhead is not distorted in graph

xcoeff=scale*xd;
ycoeff=scale*yd;

xx = [start(1), stop(1),(stop(1)+xcoeff*cos(theta+pi/4)),NaN,stop(1),... 
(stop(1)+xcoeff*cos(theta-pi/4))]';
yy = [start(2), stop(2), (stop(2)+ycoeff*sin(theta+pi/4)),NaN,stop(2),... 
(stop(2)+ycoeff*sin(theta-pi/4))]';

handle = plot(xx,yy);
