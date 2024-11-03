#include <math.h>
#include "mystruc.h" 

double bartim(double utc, int day, int year){
  double tleap();

  double tdt, g, date, n;
  double leap_seconds;

  static double epoch = -5479.5e+0;
  static double off = 32.184000e+0;
  static double sec_in_day = 86400.0e+0;
  static double t1 = 0.001658e+0;
  static double t2 = 0.000014e+0;
  static double a1 = 6.2401e+0;
  static double a2 = 0.017202e+0;

  leap_seconds = tleap(day,year);
  tdt = utc + off + leap_seconds;

  n = day + 365*(year-1985) + (year-1984)/4;
  date = n;
  date = date + epoch + (tdt/sec_in_day);

  g = a1 + a2*date;
  return(tdt + t1*sin(g) + t2*sin(2.0*g));
}
