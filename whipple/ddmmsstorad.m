function y = ddmmsstorad(x)
%
dd = fix(x/10000);
mm = fix((x - dd*10000)/100);
ss = x - dd*10000 - mm*100;
y = dd*pi/180 + mm/60*pi/180 + ss/60/60*pi/180;
