function y = upperlimit(con,coff,csigma)
%
% Usage: upperlimit(on,off,sigma)
%        on is total on counts
%        off is total off counts
%            OFF cnts for ON/OFF
%            OFF*alpha for tracking
%        sigma is the uncertainty in the excess 
%        e.g. sigma = sqrt(ON + OFF)  for ON/OFF data - Poisson error
%                   = sqrt(alpha*(ON+OFF)) for tracking with no uncertainty
%                                          in alpha ratio (OFF not scaled) -
%                                          Li and Ma
%                   = sqrt(ON + alpha^2*OFF + OFF^2*delta_alpha^2)
%                                          for tracking with uncertainty in
%                                          alpha ratio - Poisson errors and
%                                          propagation of error in alpha
disp('Helene upper limit 99.9% confidence level')
y = fmin('upperfun',2*csigma,4*csigma,[0,1.0e-20],con,coff,csigma);
disp(sprintf('   function minimum %e',upperfun(y,con,coff,csigma)));
