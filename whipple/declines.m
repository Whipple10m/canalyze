function y = declines(x)
% x and y in degrees
global ra_center dec_center ra dec
x = x*pi/180;
a = 1 + 1./x.^2/cos(dec_center)^2;
b = tan(dec_center)*tan(dec);
c = tan(dec_center)^2*tan(dec)^2 - 1./x.^2/cos(dec_center)^2;
cos_h = (-b + sqrt(b^2-4*a.*c))/2./a;
y = (cos(dec_center)*sin(dec) - sin(dec_center)*cos(dec).*cos_h)./...
    (sin(dec_center)*sin(dec) + cos(dec_center)*cos(dec).*cos_h);
y = y*180/pi;
