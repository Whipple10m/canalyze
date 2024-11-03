function y = radeclines(rac,decc,xmin,xmax,ymin,ymax,offset,color)
global ra_center dec_center ra dec
%
% ra and dec in radians
%
ra_center = rac;
dec_center = decc;
hold on
%
% declination lines
%
decdeg = fix(decc*180/pi);
decmin = fix((decc - decdeg*pi/180)*180*60/pi/30)*30;
for min = decmin+ceil(ymin*60/30)*30:30:decmin+ceil(ymax*60/30)*30
   dec = decdeg*pi/180 + min/60*pi/180;
   [x,y] = fplot('declines',[xmin xmax ymin ymax]);
   plot(x,y,color);
   axis([xmin xmax ymin ymax]);
   dec=dec+0.0001;
   str = sprintf('%2d %2d',fix(dec*180/pi),...
                  abs(round((dec*180/pi-fix(dec*180/pi))*60)));
   text(xmin-offset,(dec-decc)*180/pi,str,'FontSize',10);
end
%
% right ascension
%
rachr = fix(rac*12/pi);
racmin = fix((rac - rachr*pi/12)*12*60/pi);
for min = racmin+fix(xmin/0.25):2:racmin+fix(xmax/0.25)
   ra = rachr*pi/12 + min*pi/12/60;
   [x,y] = fplot('ralines',[xmin xmax ymin ymax]);
   plot(y,x,color);
   axis([xmin xmax ymin ymax]);
   str = sprintf('%2d %2d',fix(ra*12/pi),...
                  round((ra - fix(ra*12/pi)*pi/12)*12*60/pi));
   text((ra_center-ra)*180/pi,ymin-offset,str,'Rotation',90,'FontSize',10);
end
hold off
