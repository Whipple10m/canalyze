/****************************************************************************/
/* find_orbit: extracts binary ephemeris from data base                     */
/* Written by T.A. Hall R.W. Lessard 991109                                 */
/* Usage: find_orbit binary_name(B1950) date(yymmdd)                        */
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
void read_binary_catalog(char *, int);

int main(int argc, char *argv[])
{
  progname = (char *)basename(*argv);
  if(argc == 1){
    printf("Usage: %s binary_name(B1950) date(yymmdd)\n",progname);
    exit(EXIT_FAILURE);
  }

  read_binary_catalog(argv[1], atoi(argv[2]));

  exit(EXIT_SUCCESS);
}

void read_binary_catalog(char *source, int idate) {

  FILE *fp;
  int i;
 
  int iy, im, id;
  
  char line[160];
  char name[10];
  double mjd_begin, mjd_end, mjd;
  double epoch, Porb, Porb_dot, ecc, asini, Omega, Omega_dot, gamma;

  iy = (int)(idate/10000);
  im = (int)((idate - iy*10000)/100);
  id = (int)(idate - (iy*10000 + im*100));
  sla_caldj_(&iy,&im,&id,&mjd,&i);

  if((fp = fopen(BINARY_CATALOG,"r"))==NULL){
    printf("%s: can\'t find %s\n",progname,BINARY_CATALOG);
    exit(EXIT_FAILURE);
  }

  for(i=0;i < 5;++i){  
    fgets(line,160,fp);
  }
  
  while(fgets(line,160,fp) != NULL) {
    i = sscanf(line,"%s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
	   name,&mjd_begin,&mjd_end,&epoch,&Porb,&Porb_dot,&ecc,&asini,&Omega,
	   &Omega_dot, &gamma);
    if((!strcmp(name, source)) && (i == 11) &&
       ((mjd >= (double)mjd_begin) && (mjd <= (double)mjd_end))){
      printf("%.11lf,%.11g,%.11lf,%.11lf,%.11lf,%.11g,%.11lf,%.11lf\n",
	     Porb,Porb_dot,ecc,asini,Omega,Omega_dot,gamma,epoch);
      fclose(fp);
      exit(EXIT_SUCCESS);
    }
  }

  fclose(fp);
  return; 
}
