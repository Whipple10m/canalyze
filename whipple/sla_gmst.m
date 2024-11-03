function GMST = sla_gmst(UT1)
%   Converted to MATLAB by
%   R.W. Lessard
%   Purdue University
%   99/09/15
%      DOUBLE PRECISION FUNCTION sla_GMST (UT1)
%*+
%*     - - - - -
%*      G M S T
%*     - - - - -
%*
%*  Conversion from universal time to sidereal time (double precision)
%*
%*  Given:
%*    UT1    dp     universal time (strictly UT1) expressed as
%*                  modified Julian Date (JD-2400000.5)
%*
%*  The result is the Greenwich mean sidereal time (double
%*  precision, radians).
%*
%*  The IAU 1982 expression (see page S15 of 1984 Astronomical
%*  Almanac) is used, but rearranged to reduce rounding errors.
%*  This expression is always described as giving the GMST at
%*  0 hours UT.  In fact, it gives the difference between the
%*  GMST and the UT, which happens to equal the GMST (modulo
%*  24 hours) at 0 hours UT each day.  In this routine, the
%*  entire UT is used directly as the argument for the
%*  standard formula, and the fractional part of the UT is
%*  added separately;  note that the factor 1.0027379... does
%*  not appear.
%*
%*  See also the routine sla_GMSTA, which delivers better numerical
%*  precision by accepting the UT date and time as separate arguments.
%*
%*  Called:  sla_DRANRM
%*
%*  P.T.Wallace   Starlink   14 September 1995
%*
%*  Copyright (C) 1995 Rutherford Appleton Laboratory
%*-
D2PI=6.283185307179586476925286766559;
S2R=7.272205216643039903848711535369e-5;
%*  Julian centuries from fundamental epoch J2000 to this UT
TU=(UT1-51544.5)/36525;
%*  GMST at this UT
GMST=sla_dranrm(mod(UT1,1)*D2PI+(24110.54841+(8640184.812866+...
                (0.093104-6.2e-6*TU)*TU)*TU)*S2R);
