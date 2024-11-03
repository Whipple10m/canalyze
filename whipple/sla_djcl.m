function [IY,IM,ID,FD,J] = sla_djcl(DJM)
%   Converted to MATLAB by
%   R.W. Lessard
%   Purdue University
%   99/09/15
%      SUBROUTINE sla_DJCL (DJM, IY, IM, ID, FD, J)
%*+
%*     - - - - -
%*      D J C L
%*     - - - - -
%*
%*  Modified Julian Date to Gregorian year, month, day,
%*  and fraction of a day.
%*
%*  Given:
%*     DJM      dp     modified Julian Date (JD-2400000.5)
%*
%*  Returned:
%*     IY       int    year
%*     IM       int    month
%*     ID       int    day
%*     FD       dp     fraction of day
%*     J        int    status:
%*                       0 = OK
%*                      -1 = unacceptable date (before 4701BC March 1)
%*
%*  The algorithm is derived from that of Hatcher 1984
%*  (QJRAS 25, 53-55).
%*
%*  P.T.Wallace   Starlink   27 April 1998
%*
%*  Copyright (C) 1998 Rutherford Appleton Laboratory
%*-
%*  Check if date is acceptable
if DJM <= -2395520 | DJM >= 1e9
   J=-1;
else
   J=0;
%*     Separate day and fraction
   F=mod(DJM,1);
   if F < 0 
      F=F+1;
   end
   D=round(DJM-F);
%*     Express day in Gregorian calendar
   JD=D+2400001;
%
   N4=4*fix(JD+((6*fix((4*JD-17918)/146097))/4+1)/2-37);
   ND10=fix(10*fix(mod(N4-237,1461)/4)+5);
%
   IY=fix(N4/1461-4712);
   IM=fix(mod(ND10/306+2,12)+1);
   ID=fix(mod(ND10,306)/10+1);
   FD=F;
%
   J=0;
end
