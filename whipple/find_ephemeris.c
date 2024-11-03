/****************************************************************************/
/* find_ephemeris: extracts pulsar ephemeris from data base                 */
/* Written by R.W. Lessard 981228                                           */
/* Usage: find_ephemeris pulsar_name(B1950) date(yymmdd)                    */
/* Patchlist:                                                               */
/* Version: 1.0                                                             */
/****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifdef SUN
#include <libgen.h>
#else
#include <unistd.h>
#include <getopt.h>
#endif
#include "sla.h"

char *progname;

/*************************************************/
/* Function declarations                         */
/*************************************************/
void read_pulsar_catalog(char *, int);
void read_crab_catalog(int);

int main(int argc, char *argv[])
{
  progname = (char *)basename(*argv);
  if(argc == 1){
    printf("Usage: %s pulsar_name(B1950) date(yymmdd)\n",progname);
    exit(EXIT_FAILURE);
  }
  if(!strcmp(argv[1],"0531+21") || !strcmp(argv[1],"0534+2200"))
    read_crab_catalog(atoi(argv[2]));
  else
    read_pulsar_catalog(argv[1], atoi(argv[2]));

  exit(EXIT_FAILURE);
}

void read_pulsar_catalog(char *source, int idate) {

  FILE *fp;
  int i;
 
  float ra_s, dec_s, rms;
  int iy, im, id;
  int ra_h, ra_min, dec_deg, dec_min, dec_sgn, diff;
  
  char obs, line[160];
  char name1[10], binary, dummy[10], name2[10];
  double f, f1, f2;
  double jd0, ra, dec, mjd;
  long  jd1, jd2, mjd1, mjd2; 
  
  iy = (int)(idate/10000);
  im = (int)((idate - iy*10000)/100);
  id = (int)(idate - (iy*10000 + im*100));
  sla_caldj_(&iy,&im,&id,&mjd,&i);

  if((fp = fopen(PULSAR_CATALOG,"r"))==NULL){
    printf("%s: can\'t find %s\n",progname,PULSAR_CATALOG);
    exit(EXIT_FAILURE);
  }

  for(i=0;i < 2;++i){  
    fgets(line,160,fp);
  }
  
  while(fgets(line,160,fp) != NULL) {
    if(sscanf(line,"%s %d %d %f %d %d %f %ld %ld %lf %lf %lf %lf %f %c %c %s", 
	      name1, &ra_h, &ra_min, &ra_s, &dec_deg, &dec_min, &dec_s,
	      &jd1, &jd2, &jd0, &f, &f1, &f2, &rms, &obs, &binary, name2) == 17){
      dec_sgn = (dec_deg < 0) ? -1 : 1;
      if (binary != '*') {
	dummy[0] = binary; dummy[1] = '\0';
	strcat(dummy,name2);
	strcpy(name2,dummy);
      }
  
      if((!strcmp(name1, source) || !strcmp(name2, source)) &&
	 ((mjd >= (double)jd1) && (mjd <= (double)jd2))) {
	ra = (ra_h + ra_min/60. + ra_s/3600.) / 24 * 2 * M_PI;
	dec = dec_sgn * (abs(dec_deg) + dec_min/60. +dec_s/3600.) * M_PI / 180;
	printf("%.10lf,%.10lf,%.10lf,%.10g,%.9lf\n",ra,dec,f,f1,jd0);
	fclose(fp);
	exit(EXIT_SUCCESS);
      }
    }
  }

  fclose(fp);
  return; 
}

void read_crab_catalog(int idate)
{
  FILE *fp;
  int iy, im, id, j;
  double djm;    
  double ra,dec,freq,fdot,epoch,radio,repoch;

  iy = (int)(idate/10000);
  im = (int)((idate - iy*10000)/100);
  id = (int)(idate - (iy*10000 + im*100));
  sla_caldj_(&iy,&im,&id,&djm,&j);

  if ((fp = fopen(CRAB_CATALOG,"r")) != NULL){
    while (fscanf(fp," %lf %lf %lf %lf %lf %lf\n", 
		  &ra,&dec,&epoch,&radio,&freq,&fdot) != EOF){
      if (abs(epoch - djm) < 16.0) {
	repoch = epoch + radio/86400.0L; 
	printf("%.10lf,%.10lf,%.10lf,%.10g,%.9lf\n", 
	       ra*M_PI/12.0,dec*M_PI/180.0,freq,fdot,repoch);
	fclose(fp);
	exit(EXIT_SUCCESS);
      }
    }
  }
  fclose(fp);

  return; 
}
