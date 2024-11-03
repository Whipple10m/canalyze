function y = gupper(con,coff,alpharatio)
%
% Usage: gupper(on,off,alpharatio)
%        on is total on counts
%        off is total off counts (not scaled by alpharatio)
%        alpharatio is ratio of ON/OFF exposures (time or bins or whatever)
%
csigma=sqrt(alpharatio*(con+coff));
coff = alpharatio*coff;
y = fmin('upperfun',2*csigma,4*csigma,[0,1.0e-20],con,coff,csigma);

