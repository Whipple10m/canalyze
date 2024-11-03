function y = circ(cx,cy,rad,varargin)
%
theta = 0:0.05:6.5;
x = cos(theta);
y = sin(theta);
%
x = rad*x + cx;
y = rad*y + cy;
h=line(x,y);
set(h,varargin{:});