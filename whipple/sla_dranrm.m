function DRANRM = sla_DRANRM(ANGLE)
%   Converted to MATLAB by
%   R.W. Lessard
%   Purdue University
%   99/09/15
%      DOUBLE PRECISION FUNCTION sla_DRANRM (ANGLE)
%*+
%*     - - - - - - -
%*      D R A N R M
%*     - - - - - - -
%*
%*  Normalize angle into range 0-2 pi  (double precision)
%*
%*  Given:
%*     ANGLE     dp      the angle in radians
%*
%*  The result is ANGLE expressed in the range 0-2 pi (double
%*  precision).
%*
%*  P.T.Wallace   Starlink   23 November 1995
%*
%*  Copyright (C) 1995 Rutherford Appleton Laboratory
%*-
D2PI=6.283185307179586476925286766559;
DRANRM=mod(ANGLE,D2PI);
if DRANRM < 0
   DRANRM = DRANRM + D2PI;
end
