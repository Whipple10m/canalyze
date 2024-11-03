function h = circle(cx,cy,radius,color)
% FUNCTION H = CIRCLE(X,Y,R,C)
% Plots a circle with radius R centerd at X,Y with fill color 
% specified by RGB color C
% Returns handle to fill
% Written by R.W. Lessard
%            Purdue University
%            03/09/99
%
t = 0:.1:2*pi;
x = radius*cos(t) + cx;
y = radius*sin(t) + cy;
h = fill(x,y,color);
