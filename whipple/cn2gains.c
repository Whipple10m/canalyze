/***************************************************************************/
/* cn2gains: calculates nitrogen gains from hdf files                      */
/* Written by R.W. Lessard 961128                                          */
/* Usage: cn2gains [options] *.hdf                                         */
/* Patchlist:                                                              */
/* Version: HDF                                                            */
/***************************************************************************/
/*************************************************************************/
/* includes                                                              */
/*************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef LINUX
#include <unistd.h>
#include <getopt.h>
#endif
#include <libgen.h>
#include <math.h>
#include "hdf.h"

/*************************************************************************/
/* defines                                                               */
/*************************************************************************/
#define EXPAND(argc, argv)
#define DEFAULT_PATH "./"
#define FOREVER         1
#define MINADC 50.0
#define MAXADC 1000.0
#define MINPERCENT 0.75

/*************************************************************************/
/* global variables                                                      */
/*************************************************************************/
char *filename;
char *progname;

/*************************************************************************/
/* command line options                                                  */
/*************************************************************************/
char path_to_db[FILENAME_MAX];
extern int optind;

/*************************************************************************/
/* function declarations                                                 */
/*************************************************************************/
void getoptions();
void usage();

main(int argc, char **argv)
{
  char fileid[9], outname[FILENAME_MAX];
  double utc_start, dfr;
  FILE *outfile;
  float *peds, *pedvars, *gains, *fadc, meanadc;
  int *ngains, nadc, npmt, version, nmean, nhigh, nlow;
  int j, imm, idd, iyy, idate;
  short code, *adc;

  int32 file_id, vevent_ref, vevent_id, vruninfo_ref, vruninfo_id;
  uchar8 *runinfobuf, *eventbuf;

  EXPAND(argc, argv);
  progname = basename(*argv); progname++;
  if(argc == 1){
    usage();
    exit(EXIT_FAILURE);
  }

/*************************************************************************/
/* set default command line options                                      */
/*************************************************************************/
  strcpy(path_to_db,DEFAULT_PATH);

/*************************************************************************/
/* get command line options                                              */
/*************************************************************************/
  getoptions(argc, argv);

/*************************************************************************/
/* MAIN LOOP                                                             */
/*************************************************************************/
  while(optind < argc){
    filename = *(argv+optind);
    ++optind;

/*************************************************************************/
/* open HDF file                                                         */
/*************************************************************************/
    file_id = Hopen(filename,DFACC_RDONLY,0);
    if(file_id == -1){
      fprintf(stderr,"%s: can\'t open %s\n",progname,filename);
      exit(EXIT_FAILURE);
    }
    Vstart(file_id);

/*************************************************************************/
/* set runinfo fields to read                                            */
/*************************************************************************/
    vruninfo_ref = VSfind(file_id,"10M Run Information");
    vruninfo_id = VSattach(file_id,vruninfo_ref,"r");
    VSsetfields(vruninfo_id,"VERSION,NADC,UTC_START");
    runinfobuf = (uchar8 *)malloc(2*sizeof(int32)+sizeof(float64));

/*************************************************************************/
/* read run information                                                  */
/*************************************************************************/
    VSread(vruninfo_id,runinfobuf,1,FULL_INTERLACE);
    memcpy(&version,runinfobuf,4);
    memcpy(&nadc,runinfobuf+4,4);
    memcpy(&utc_start,runinfobuf+8,8);

/*************************************************************************/
/* calculate run date and output run information                         */
/*************************************************************************/
    strncpy(fileid,basename(filename),8); fileid[8] = '\0';
    printf("%s: reading %s ...\n",progname,filename);
    sla_djcl_(&utc_start,&iyy,&imm,&idd,&dfr,&j);
    if(iyy < 2000)
      idate = (iyy-1900)*10000+imm*100+idd;
    else
      idate = (iyy-2000)*10000+imm*100+idd;
    printf("   run date: %06d\n",idate);
    printf("   no. channels: %d\n",nadc);
    npmt = get_npmt(nadc);
    printf("   no. pmts:     %d\n",npmt);

/*************************************************************************/
/* set event fields to read                                              */
/*************************************************************************/
    vevent_ref = VSfind(file_id,"10M Event");
    vevent_id = VSattach(file_id,vevent_ref,"r");
    VSsetfields(vevent_id,"CODE,ADC");
    eventbuf = (uchar8 *)malloc(sizeof(int16) + nadc*sizeof(int16));

/*************************************************************************/
/* allocate memory and initialize arrays                                 */
/*************************************************************************/
    adc = (short *)malloc(nadc*sizeof(short));
    fadc = (float *)malloc(nadc*sizeof(float));
    peds = (float *)malloc(nadc*sizeof(float));
    pedvars = (float *)malloc(nadc*sizeof(float));
    gains = (float *)malloc(nadc*sizeof(float));
    ngains = (int *)malloc(nadc*sizeof(int));
    for(j = 0;j < nadc;j++){
      gains[j] = 0.0;
      ngains[j] = 0;
    }
/*************************************************************************/
/* see if it's been done                                                 */
/*************************************************************************/
    if(getn2gains(fileid,idate,gains,path_to_db)){
      printf("   n2gains already exist, skipping file ...\n");
    }else if(!getpeds(fileid,idate,peds,pedvars,path_to_db)){
      printf("   can\'t find pedestals, skipping file ...\n");
    }else{
/*************************************************************************/
/* read until the end of the file                                        */
/*************************************************************************/
      while(VSread(vevent_id,eventbuf,1,FULL_INTERLACE)==1){
	memcpy(&code,eventbuf,2);
	if(code & 4){
	  memcpy(adc,eventbuf+2,2*nadc);
	  meanadc = 0.0; nmean = 0; nlow = 0; nhigh = 0;
	  for(j = 0;j < npmt;j++){
	    fadc[j] = (float)adc[j] - peds[j];
	    if(fadc[j] < MINADC)
	      ++nlow;
	    else if(fadc[j] > MAXADC)
	      ++nhigh;
	    else{
	      meanadc += fadc[j];
	      nmean++;
	    }
	  }
	  if(nmean > (int)(npmt*MINPERCENT)){
	    meanadc /= (float)nmean;
	    for(j = 0;j < npmt;j++){
	      if((fadc[j] > MINADC) && (fadc[j] < MAXADC)){
		gains[j] += meanadc/fadc[j];
		ngains[j]++;
	      }
	    }
	  }else{
	    printf("   ** warning ** %3d saturated %3d too low\n",
		   nhigh,nlow);
	  }
	}
      }
/*************************************************************************/
/* calculate gains                                                       */
/*************************************************************************/
      for(j = 0;j < nadc;j++){
	if(ngains[j] > 0){
	  gains[j] /= ngains[j];
	}else{
	  gains[j] = 1.0;
	  if(j < npmt){
	    printf("   ** warning ** pmt %3d has insufficient data\n",j+1);
	  }
	}
      }
/*************************************************************************/
/* output results                                                        */
/*************************************************************************/
      sprintf(outname,"%shrc%02d.cn2gains",path_to_db,(int)(idate/10000));
      if((outfile = fopen(outname, "a")) == NULL){
	printf("   **error** cannot open %s\n",outname);
	exit(EXIT_FAILURE);
      }
      fprintf(outfile,"%8s %06d %6d\n",fileid,idate,nadc);
      for(j = 0; j < nadc;j++){
	fprintf(outfile,"%7.3f",gains[j]);
	if(((j+1)%10) == 0)
	  fprintf(outfile,"\n");
      }
      if((j%10) != 0)
	fprintf(outfile,"\n");
      fclose(outfile);
    }
/*************************************************************************/
/* close HDF file and free buffers                                       */
/*************************************************************************/
    VSdetach(vruninfo_id);
    VSdetach(vevent_id);
    Vend(file_id);
    Hclose(file_id);
    free(runinfobuf); free(eventbuf); 
    free(fadc); free(adc); 
    free(peds); free(pedvars); free(gains); free(ngains);
  }
/*************************************************************************/
/* next file                                                             */
/*************************************************************************/
  printf("done ...\n");
}

/*************************************************************************/
/* function: getopt                                                      */
/* passed:   number of command line arguments  (int)                     */
/*           pointer to command line arguments (char**)                  */
/* returns:  nothing                                                     */
/* purpose:  gets command line options                                   */
/*************************************************************************/
void getoptions(int argc, char **argv)
{
  extern char *optarg;
  int c;

  while(FOREVER){
/*************************************************/
/* Note: getopt is POSIX.1 compliant thus should */
/*       available on all newer workstations.    */
/*       getopt_long would be better as it sorts */
/*       argv into options -> non-options, hence */
/*       the options could appear anywhere on    */
/*       the command line                        */
/*************************************************/  
    c = getopt(argc, argv, "g:h");
    switch(c){
    case 'g':
      strcpy(path_to_db,optarg);
      break;
    case 'h':
      usage();
      exit(EXIT_SUCCESS);
    case -1:
      return;
    default:
      printf("%s: ignoring unknown option %c\n",progname,c);
      break;
    }
  }
}

/*************************************************************************/
/* function: usage                                                       */
/* passed:   nothing                                                     */
/* returns:  nothing                                                     */
/* purpose:  prints usage information                                    */
/*************************************************************************/
void usage(void)
{
  printf("%s [options] filename.hdf\n",progname);
  printf(" where [options] are\n");
  printf("   -g %%s  path to database files [%s]\n",DEFAULT_PATH);
  printf("   -h     help\n");

  return;
}
/* Thats all folks */
