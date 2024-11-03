function [newra,newdec] = precessfrom2000(ra,dec,to)
%
T = (to - 2000)/100;
M = 1.2812323*pi/180*T + 0.0003879*pi/180*T^2 + 0.0000101*pi/180*T^3;
N = 0.5567530*pi/180*T - 0.0001185*pi/180*T^2 - 0.0000116*pi/180*T^3;
%
ram = ra - 1/2*(M + N*sin(ra).*tan(dec));
decm = dec - 1/2*N*cos(ram);
%
newra = ra + M + N*sin(ram).*tan(decm);
newdec = dec + N*cos(ram);
