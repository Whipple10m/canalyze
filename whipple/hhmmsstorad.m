function y = hhmmsstorad(x)
%
hh = fix(x/10000);
mm = fix((x - hh*10000)/100);
ss = x - hh*10000 - mm*100;
y = hh*pi/12 + mm*pi/12/60 + ss*pi/12/60/60;
