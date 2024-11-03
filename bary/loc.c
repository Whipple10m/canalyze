#include "mystruc.h"

loc(lwrtec,lcoords,gmst)

VECTOR *lwrtec;
VECTOR *lcoords;
double gmst;

{
	double st,xl,zl,rad,rap;

	double pow();
	double sqrt();
	double sin();
	double cos();

	static double conv = 7.2722052e-5;
	static double radian = 57.29578e+0;
	static double rec = 0.02127653e+0;
	static double c = 2.9979246e+8;
	static double flat = 0.00335281e+0;

	st = conv*(3600.0e+0*(lcoords->xcomp) + gmst);
	xl = (lcoords->ycomp)/radian;

	zl = sin(xl);
	rad = rec/sqrt(1.0e+0-((zl*zl)*(1.0e+0-
                       ((1.0-flat)*(1.0-flat)))));
	rap = rad * ((1.0-flat)*(1.0-flat));
	rad = rad + (lcoords->zcomp)/c;
	rap = rap + (lcoords->zcomp)/c;

	lwrtec->zcomp = rap * zl;
	zl = cos(xl);
	lwrtec->ycomp = rad * zl * sin(st);
	lwrtec->xcomp = rad * zl * cos(st);

}






