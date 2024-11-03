function y = ellipse(length,width,slope,xc,yc)
ecc = sqrt(1.0-(width/length)^2); 
slr = width^2/length; 
delta =atan(slope);
polex = xc+(length-(slr/(1.0+ecc)))*cos(delta); 
poley = yc+(length-(slr/(1.0+ecc)))*sin(delta); 
%
theta = 0.; 
for i=1:72
   if theta <= 360.
      r = slr/(1.0+ecc*cos(theta*pi/180.0-delta)); 
      xcoord(i) = r*cos(theta*pi/180.0)+polex; 
      ycoord(i) = r*sin(theta*pi/180.0)+poley; 
   end 
   theta = theta + 5.0; 
end
plot(xcoord(1:72),ycoord(1:72),'k-')
 








