function y = radtohhmmss(x)
%
hr = x*12/pi;
min = (hr - fix(hr))*60;
sec = (min - fix(min))*60;
y = fix(hr)*10000+fix(min)*100+sec;
