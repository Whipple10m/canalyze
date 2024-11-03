function y = upperfun(x,on,off,sigma)
%
y = abs(erfc(((x-(on-off))/sigma)/sqrt(2))/2 - ...
    0.001*erfc((-(on-off)/sigma)/sqrt(2))/2);
