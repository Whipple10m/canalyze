#include "mystruc.h"
           
unit(unit_vector,ra,dec)

double ra,dec;
VECTOR *unit_vector;

{
	double a,d;

	double sin();
	double cos();

	static double twopi = 6.283185308e+0;
	static double radian = 57.2957795e+0;

	a = twopi * ra/24.0e+0;
	d = dec/radian;

	unit_vector->xcomp = cos(d);
	unit_vector->ycomp = unit_vector->xcomp;
	unit_vector->xcomp = unit_vector->xcomp * cos(a);
	unit_vector->ycomp = unit_vector->ycomp * sin(a);
	unit_vector->zcomp = sin(d);


}
