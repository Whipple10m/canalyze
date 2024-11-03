#include "mystruc.h"

void earth(VECTOR *earth_position, double tdb, int day, int year, 
	   POSITION table[])
{
  int i,day_step;
  double c[4];
  double factor;
  static double au = 499.0047837L;
  static double sec_in_day = 86400.0L;

  factor = tdb/sec_in_day;

  c[0] = -1.0*factor*(factor-1.0L)*(factor-2.0L)/6.0L;
  c[1] = 0.5L*((factor*factor)-1.0L)*(factor-2.0L);
  c[2] = -0.5L*factor*(factor+1.0L)*(factor-2.0L);
  c[3] = factor*((factor*factor)-1.0L)/6.0L;

  earth_position->xcomp = 0.0L;
  earth_position->ycomp = 0.0L;
  earth_position->zcomp = 0.0L;

  day_step = day-1;

  for(i=0; i<=3; i++){
    earth_position->xcomp = earth_position->xcomp +
      c[i]*table[day_step].x;
    
    earth_position->ycomp = earth_position->ycomp +
      c[i]*table[day_step].y;
	  
    earth_position->zcomp = earth_position->zcomp +
      c[i]*table[day_step].z;
	  
    day_step = day_step + 1;
  }

  earth_position->xcomp = au*(earth_position->xcomp);
  earth_position->ycomp = au*(earth_position->ycomp);
  earth_position->zcomp = au*(earth_position->zcomp);

  return;

}
