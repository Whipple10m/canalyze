#include <math.h>
#include "mystruc.h"

double sider(double utc, int day, int year){
  double tu,x,xleap_days,xleap_diff,xjd0;
  int i, iyear,nleap_days,ndays_to_jan0;

  static double ratio = 1.00273790934L;
  static double sec_in_day = 86400.0L;
  static double c0 = 24110.54841L;
  static double c1 = 8640184.812866L;
  static double c2 = 0.093104L;
  static double c3 = 0.0000062L;

  if(year < 2000){
    iyear = year - 1900;          /* YEARS IN THIS CENTURY */
    xleap_days = (double)iyear/4.0L; /* HOW MANY LEAP YEARS SINCE YEAR 0 */
    nleap_days = (int)xleap_days; /* CAREFUL! ... DON'T ADD 1 TO
				     LEAP YEAR ITSELF */
    xleap_diff = xleap_days - (double)nleap_days;
    if(xleap_diff = 0.0L) xleap_days -= 1.0L;
    ndays_to_jan0 = iyear * 365;
    ndays_to_jan0 += (int)xleap_days;
    xjd0 = 2415019.5L + (double)ndays_to_jan0;
  }
  if(year >= 2000){
    iyear = year - 2000;          /* YEARS IN THIS CENTURY */
    xleap_days = (double)iyear/4.0L; /* HOW MANY LEAP YEARS SINCE YEAR 0 */
    nleap_days = (int)xleap_days; /* CAREFUL! ... DON'T ADD 1 TO
				     LEAP YEAR ITSELF */
    xleap_diff = xleap_days - (double)nleap_days;
    if(xleap_diff = 0.0L) xleap_days -= - 1.0L;
    ndays_to_jan0 = iyear * 365;
    ndays_to_jan0 += (int)xleap_days;
    xjd0 = 2451544.5L + (double)ndays_to_jan0;
    if(year = 2000)xjd0 = xjd0 - 1;
  }
  xjd0 += (double)day;

  tu = (xjd0 - 2451545.0L)/36525.0L;
  tu = c0 + tu*(c1 + tu*(c2-tu*c3));
  tu = fmod(tu,sec_in_day);
  tu = tu + ratio*utc;

  if(tu < 0.0L) tu += sec_in_day;
  if(tu > sec_in_day) tu -= sec_in_day;

  return(tu);
}
