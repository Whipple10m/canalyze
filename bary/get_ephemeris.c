#include "mystruc.h"

void get_ephemeris(int year, POSITION table[])
{
  FILE* in;
  int i,n,m;
  double jd, jd84;
  double day,hx,hy,hz,hvx,hvy,hvz;

  jd84 = 2445699.5L;

  n = year-1984;
  m = (n+3)/4;
  jd = jd84 + 365*n + m;

  if((in = fopen("/home/lessard/bin/earth.dat","r")) == NULL){
    printf("could not open ephemeris file\n");
    exit(-1);
  }

  jd = jd - 1.01;
  fscanf(in,"%lf %lf %lf %lf %lf %lf %lf\n",
	 &day,&hx,&hy,&hz,&hvx,&hvy,&hvz);

  while(day < jd){
    fscanf(in,"%lf %lf %lf %lf %lf %lf %lf\n",
	   &day,&hx,&hy,&hz,&hvx,&hvy,&hvz);
  }

  for(i=0; i<=ROWMAX-1; i++){
    fscanf(in,"%lf %lf %lf %lf %lf %lf %lf\n",
	   &day,
	   &table[i].x,
	   &table[i].y,
	   &table[i].z,
	   &table[i].vx,
	   &table[i].vy,
	   &table[i].vz);
  }
  fclose(in);

  return;
}
