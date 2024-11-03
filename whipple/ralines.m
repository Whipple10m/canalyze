function y = ralines(x)
% x and y in degrees
global ra_center dec_center ra dec
x = dec_center + x*pi/180;
y = (cos(x)*sin(ra_center - ra))./...
    (sin(dec_center)*sin(x) + cos(dec_center)*cos(x)*cos(ra_center - ra));
y = y*180/pi;