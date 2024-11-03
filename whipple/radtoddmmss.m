function y = radtoddmmss(x)
%
deg = x*180/pi;
min = (deg - fix(deg))*60;
sec = (min - fix(min))*60;
y = fix(deg)*10000+fix(min)*100+sec;
