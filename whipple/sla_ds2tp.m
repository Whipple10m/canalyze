function [XI, ETA, J] = sla_ds2tp(RA, DEC, RAZ, DECZ)
%      SUBROUTINE sla_DS2TP (RA, DEC, RAZ, DECZ, XI, ETA, J)
%*+
%*     - - - - - -
%*      D S 2 T P
%*     - - - - - -
%*
%*  Projection of spherical coordinates onto tangent plane:
%*  "gnomonic" projection - "standard coordinates" (double precision)
%*
%*  Given:
%*     RA,DEC      dp   spherical coordinates of point to be projected
%*     RAZ,DECZ    dp   spherical coordinates of tangent point
%*
%*  Returned:
%*     XI,ETA      dp   rectangular coordinates on tangent plane
%*     J           int  status:   0 = OK, star on tangent plane
%*                                1 = error, star too far from axis
%*                                2 = error, antistar on tangent plane
%*                                3 = error, antistar too far from axis
%*
%*  P.T.Wallace   Starlink   18 July 1996
%*
%*  Copyright (C) 1996 Rutherford Appleton Laboratory
%*-
TINY=1e-6;
%*  Trig functions
SDECZ=sin(DECZ);
SDEC=sin(DEC);
CDECZ=cos(DECZ);
CDEC=cos(DEC);
RADIF=RA-RAZ;
SRADIF=sin(RADIF);
CRADIF=cos(RADIF);
%*  Reciprocal of star vector length to tangent plane
DENOM=SDEC.*SDECZ+CDEC.*CDECZ.*CRADIF;
%*  Handle vectors too far from axis
J=zeros(1,length(RA))+3;
J(find(DENOM>TINY))=0;
J(find(DENOM>=0 & DENOM<=TINY))=1;
DENOM(find(DENOM>=0 & DENOM<=TINY))=TINY;
J(find(DENOM>-TINY & DENOM>=0 & DENOM<=TINY))=2;
DENOM(find(DENOM>-TINY & DENOM>=0 & DENOM<=TINY))=-TINY;
%*  Compute tangent plane coordinates (even in dubious cases)
XI=CDEC.*SRADIF./DENOM;
ETA=(SDEC.*CDECZ-CDEC.*SDECZ.*CRADIF)./DENOM;
