#include "mystruc.h"

/*********************************************************************/
/* function bary                                                     */
/* utc - universal time coordinated since 00hr in seconds            */
/* day - day number                                                  */
/* year - year in the form yyyy                                      */
/* right_ascension - source right ascension in decimal hours         */
/* declination - source declination in decimal degrees               */
/*********************************************************************/
double bary(double utc, int day, int year, double right_ascension,
	    double declination){

  /* table consisting of x,y,z vx,vy,vz of earth at 0hr */
  static POSITION table[ROWMAX]; 

  /* geodetic position of lab and altitude */
  /* lab_coords[1] -> west long in hours */
  /* lab_coords[2] -> latitude in degrees */
  /* altitude in meters above mean sea level */

  VECTOR lab_coords;

  /* lab position with respect to the earth's center */
  
  VECTOR lab_wrt_earth_center;

  /* earths position in rectangular celestial coord. */

  VECTOR earth_position;

  /* unit vector in source direction */
  
  VECTOR unit_vector;

  double dot_product;		/* dot product */
  double gmst;			/* Greenwich Mean Time */

  static double seconds_in_day = 86400.0L;	/* seconds in a day */
  static int first_time=1;

  double tdb;		/* Solar Barycentric Dynamic Time */

  double bartim();
  double sider();

  lab_coords.xcomp = -7.3922;
  lab_coords.ycomp = 31.683;
  lab_coords.zcomp = 2336.0;

  if(first_time){
    get_ephemeris(year,table);
    first_time=0;
  }

  tdb = bartim(utc,day,year);
  earth(&earth_position,tdb,day,year,table);
  gmst = sider(utc,day,year);
  loc(&lab_wrt_earth_center,&lab_coords,gmst);
  unit(&unit_vector,right_ascension,declination);

  dot_product = unit_vector.xcomp *
    (lab_wrt_earth_center.xcomp + earth_position.xcomp);

  dot_product = dot_product + unit_vector.ycomp *
    (lab_wrt_earth_center.ycomp + earth_position.ycomp);
  
  dot_product = dot_product + unit_vector.zcomp *
    (lab_wrt_earth_center.zcomp + earth_position.zcomp);

  if (dot_product > 500.0){
    printf("   ** error ** light travel time b/w earth/SSB");
    printf("is greater than  maximum possible\n");   
  }
  return(dot_product + tdb);
}
  
