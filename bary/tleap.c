double tleap(int day, int year){
  double leap_seconds;

  if((year == 1985) && (day <= 181)) leap_seconds=22.0e+0;
  if((year == 1985) && (day > 181)) leap_seconds=23.0e+0;
  if((year == 1986) || (year == 1987)) leap_seconds=23.0e+0;
  if(year > 1987) leap_seconds=24.0e+0;
  if(year == 1990) leap_seconds=25.0e+0;
  if(year == 1991) leap_seconds=26.0e+0;
  if((year == 1992) && (day <= 181)) leap_seconds=26.0e+0;
  if((year == 1992) && (day > 181)) leap_seconds=27.0e+0;
  if((year == 1993) && (day <= 181)) leap_seconds=27.0e+0;
  if((year == 1993) && (day > 181)) leap_seconds=28.0e+0;
  if((year == 1994) && (day <= 181)) leap_seconds=28.0e+0;
  if((year == 1994) && (day > 181)) leap_seconds=29.0e+0;
  if(year == 1995) leap_seconds=29.0e+0;
  if(year == 1996) leap_seconds=30.0e+0;
  if((year == 1997) && (day <= 181)) leap_seconds=30.0e+0;
  if((year == 1997) && (day > 181)) leap_seconds=31.0e+0;
  if(year == 1998) leap_seconds=31.0e+0;
  if(year == 1999) leap_seconds=32.0e+0;
	
  return(leap_seconds);
}
