/***************************************************************************/
/* cpeds: calculates pedestals from hdf files                              */
/* Written by R.W. Lessard 961128                                          */
/* Usage: cpeds [options] *.hdf                                            */
/* Patchlist:                                                              */
/* Version: HDF                                                            */
/*          Y2K compliant 12/27/1999                                       */
/*          MPI 02/21/00                                                   */
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
#ifdef MPI
#include "mpi.h"
#endif

/*************************************************************************/
/* defines                                                               */
/*************************************************************************/
#define DEFAULT_PATH "./"
#define MAXADC 100
#define MAXPED 50.0
#define MINPED 10.0
#define MAXPEDVAR 7.0
#define MINPEDVAR 2.0

/*************************************************************************/
/* global variables                                                      */
/*************************************************************************/
char *filename;
char *progname;
int myid = 0, master = 0;

/*************************************************************************/
/* command line options                                                  */
/*************************************************************************/
char path_to_db[FILENAME_MAX];
float upper_threshold;
float lower_threshold;
extern int optind;

/*************************************************************************/
/* function declarations                                                 */
/*************************************************************************/
void getoptions();
void usage();

int main(int argc, char **argv)
{
  char fileid[9], coffname[FILENAME_MAX], cpedsname[FILENAME_MAX];
  double utc_start, dfr;
  FILE *cofffile, *cpedsfile;
  float *peds, *pedvars, *pedsbuf, *pedvarsbuf;
  int *npeds, *npedsbuf, nadc, npmt, version;
  int i, j, k, imm, idd, iyy, idate;
  short code, *adc;

  int32 file_id, vevent_ref, vevent_id, vruninfo_ref, vruninfo_id;
  uchar8 *runinfobuf, *eventbuf;

  int numproc = 1, ierr = 0, proc = 1;

#ifdef MPI
/*************************************************************************/
/* initialize MPI                                                        */
/*************************************************************************/
  MPI_Status stat;
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numproc);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
#endif

  progname = basename(*argv);
  if(argc == 1){ /* no arguments print usage and exit */
    if(myid == master)
      usage();
#ifdef MPI
    MPI_Finalize();
#endif
    return 1;
  }
/*************************************************************************/
/* set default command line options                                      */
/*************************************************************************/
  strcpy(path_to_db,DEFAULT_PATH);
  upper_threshold = MAXPEDVAR;
  lower_threshold = MINPEDVAR;

/*************************************************************************/
/* get command line options                                              */
/*************************************************************************/
  getoptions(argc, argv);
  filename = *(argv+optind);

  if(myid == master){ /* master does IO */
/*************************************************************************/
/* open HDF file                                                         */
/*************************************************************************/
    file_id = Hopen(filename,DFACC_RDONLY,0);
    if(file_id == -1){
      fprintf(stderr,"%s: can\'t open %s\n",progname,filename);
      ierr = 1; /* set error flag */
    }else{
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

      if(getpeds(fileid,idate,NULL,NULL,path_to_db)){
	printf("   pedestals already exist, skipping file ...\n");
	ierr = 1; /* set error flag */
      }else{
/*************************************************************************/
/* set event fields to read                                              */
/*************************************************************************/
	vevent_ref = VSfind(file_id,"10M Event");
	vevent_id = VSattach(file_id,vevent_ref,"r");
	VSsetfields(vevent_id,"CODE,ADC");
	eventbuf = (uchar8 *)malloc(sizeof(int16) + nadc*sizeof(int16));
      }
    }
  }

#ifdef MPI
  MPI_Bcast(&nadc,1,MPI_INT,master,MPI_COMM_WORLD); /* send nadc to all */
  MPI_Bcast(&ierr,1,MPI_INT,master,MPI_COMM_WORLD); /* send ierr to all */
#endif
  if(ierr){ /* quit on error */
#ifdef MPI
    MPI_Finalize();
#endif
    return 1;
  }

/*************************************************************************/
/* allocate memory and initialize arrays                                 */
/*************************************************************************/
  adc = (short *)malloc(nadc*sizeof(short));
  peds = (float *)malloc(sizeof(float)*nadc);
  pedsbuf = (float *)malloc(sizeof(float)*nadc);
  pedvars = (float *)malloc(sizeof(float)*nadc);
  pedvarsbuf = (float *)malloc(sizeof(float)*nadc);
  npeds = (int *)malloc(sizeof(int)*nadc);
  npedsbuf = (int *)malloc(sizeof(int)*nadc);
  for(j = 0;j < nadc;j++){
    peds[j] = 0.0;
    pedvars[j] = 0.0;
    npeds[j] = 0;
  }

/*************************************************************************/
/* MAIN LOOP                                                             */
/*************************************************************************/
  while(1){
/*************************************************************************/
/* master                                                                */
/*************************************************************************/
    if(myid == master){ /* master reads file and sends data to slaves */
      if(VSread(vevent_id,eventbuf,1,FULL_INTERLACE)==1){
	memcpy(&code,eventbuf,2);
	if(code & 1){
	  memcpy(adc,eventbuf+2,2*nadc);
#ifdef MPI	  
	  if(numproc > 1){
	    MPI_Send(adc,nadc,MPI_SHORT,proc,1,MPI_COMM_WORLD);
	    if(proc < numproc-1) proc++; else proc = 1;
	  }
#endif
	}
      }else{
#ifdef MPI	
	for(j=1;j < numproc;j++){ /* done send stop and get results */
	  MPI_Send(adc,0,MPI_SHORT,j,0,MPI_COMM_WORLD);
	  MPI_Recv(pedsbuf,nadc,MPI_FLOAT,j,1,MPI_COMM_WORLD,&stat);
	  MPI_Recv(pedvarsbuf,nadc,MPI_FLOAT,j,2,MPI_COMM_WORLD,&stat);
	  MPI_Recv(npedsbuf,nadc,MPI_INT,j,3,MPI_COMM_WORLD,&stat);
	  for(i=0;i < nadc;i++){
	    peds[i] += pedsbuf[i];
	    pedvars[i] += pedvarsbuf[i];
	    npeds[i] += npedsbuf[i];
	  }
	}
#endif
/*************************************************************************/
/* calculate pedestals and output some information                       */
/*************************************************************************/
	for(j = 0;j < nadc;j++){
	  if(npeds[j] >= 2){
	    peds[j] /= npeds[j];
	    pedvars[j] = sqrt((pedvars[j] - (float)(npeds[j])*peds[j]*peds[j])/
			      ((float)(npeds[j]) - 1.0));
	    if(((peds[j] < MINPED) || (peds[j] > MAXPED)) && (j < npmt))
	      printf("   ** info ** pmt %3d has a pedestal of %6.2f\n",
		     j+1,peds[j]);
	    if(((pedvars[j] < MINPEDVAR) || (pedvars[j] > MAXPEDVAR))
	       && (j < npmt)){
	      printf("   ** info ** pmt %3d has a pedestal variance of %6.2f\n"
		     ,j+1,pedvars[j]);
	    }
	  }else{
	    peds[j] = 0.0;
	    pedvars[j] = 0.0;
	    if(j < npmt)
	      printf("   ** warning ** pmt %3d has insufficient data\n",j+1);
	  }
	}
/*************************************************************************/
/* output results                                                        */
/*************************************************************************/
	sprintf(cpedsname,"%shrc%02d.cpeds",path_to_db,(int)(idate/10000));
	sprintf(coffname,"%shrc%02d.coff",path_to_db,(int)(idate/10000));
	if((cpedsfile = fopen(cpedsname, "a")) == NULL){
	  fprintf(stderr,"   **error** can\'t open %s\n",cpedsname);
#ifdef MPI
	  MPI_Finalize();
#endif
	  return 1;
	}
	if((cofffile = fopen(coffname, "a")) == NULL){
	  fprintf(stderr,"   **error** can\'t open %s\n",coffname);
#ifdef MPI
	  MPI_Finalize();
#endif
	  return 1;
	}

/*************************************************************************/
/* peds                                                                  */
/*************************************************************************/
	fprintf(cpedsfile,"%8s %06d %6d\n",fileid,idate,nadc);
	for(j = 0; j < nadc;j++){
	  fprintf(cpedsfile,"%7.3f",peds[j]);
	  if(((j+1)%10) == 0)
	    fprintf(cpedsfile,"\n");
	}
	if((j%10) != 0)
	  fprintf(cpedsfile,"\n");

/*************************************************************************/
/* pedvars                                                               */
/*************************************************************************/
	for(j = 0; j < nadc;j++){
	  fprintf(cpedsfile,"%7.3f",pedvars[j]);
	  if(((j+1)%10) == 0)
	    fprintf(cpedsfile,"\n");
	}
	if((j%10) != 0)
	  fprintf(cpedsfile,"\n");

/*************************************************************************/
/* tubes off                                                             */
/*************************************************************************/
	for(i = 0, k = 0;i < npmt;i++){
	  if((pedvars[i] >= upper_threshold) || 
	     (pedvars[i] <= lower_threshold)){
	    if(k == 0)
	      fprintf(cofffile,"%10s",fileid);
	    fprintf(cofffile,"%5d",i+1);
	    k++;
	  }
	}
	if(k != 0) fprintf(cofffile,"\n");

	fclose(cpedsfile); fclose(cofffile);

/*************************************************************************/
/* close HDF file and return                                             */
/*************************************************************************/
	VSdetach(vruninfo_id); VSdetach(vevent_id);
	Vend(file_id); Hclose(file_id);
#ifdef MPI
	MPI_Finalize();
#endif
	printf("done ...\n");
	
	return 0;
      }
    }
/*************************************************************************/
/* slave                                                                 */
/*************************************************************************/
    if(myid != 0 || (numproc == 1 && code & 1)){
#ifdef MPI
      if(numproc > 1){
	MPI_Recv(adc,nadc,MPI_SHORT,master,
		 MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
	if(!stat.MPI_TAG){ /* received stop, send results */
	  MPI_Send(peds,nadc,MPI_FLOAT,master,1,MPI_COMM_WORLD);
	  MPI_Send(pedvars,nadc,MPI_FLOAT,master,2,MPI_COMM_WORLD);
	  MPI_Send(npeds,nadc,MPI_INT,master,3,MPI_COMM_WORLD);
	  MPI_Finalize();
	  return 0;
	}
      }
#endif
/*************************************************************************/
/* accumulate sums for average and variance                              */
/*************************************************************************/
      for(j = 0;j < nadc;j++){
	if((abs(adc[j]) < MAXADC) && (abs(adc[j]) > 0)){
	  peds[j] += (float)(abs(adc[j]));
	  pedvars[j] += (float)(abs(adc[j]))*(float)(abs(adc[j]));
	  ++npeds[j];
	}
      }
    }
  }
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

  while(1){
/*************************************************/
/* Note: getopt is POSIX.1 compliant thus should */
/*       available on all newer workstations.    */
/*       getopt_long would be better as it sorts */
/*       argv into options -> non-options, hence */
/*       the options could appear anywhere on    */
/*       the command line                        */
/*************************************************/  
    c = getopt(argc, argv, "g:hl:u:");
    switch(c){
    case 'g':
      strcpy(path_to_db,optarg);
      break;
    case 'h':
      if(myid == master)
	usage();
#ifdef MPI
      MPI_Finalize();
#endif
      exit(EXIT_SUCCESS);
    case 'l':
      if(sscanf(optarg,"%f",&lower_threshold) != 1)
	printf("%s: invalid lower threshold option\n",progname);
      break;
    case 'u':
      if(sscanf(optarg,"%f",&upper_threshold) != 1)
	printf("%s: invalid upper threshold option\n",progname);
      break;
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
  printf("   -l %%f  lower pedestal variance threshold [%.2f]\n",
	 MINPEDVAR);
  printf("   -u %%f  upper pedestal variance threshold [%.2f]\n",
	 MAXPEDVAR);   
  return;
}
/* Thats all folks */
