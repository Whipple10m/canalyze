/****************************************************************************/
/* red2hdf: converts reduced files to HDF format                            */
/* Written by R.W. Lessard 000120                                           */
/* Usage: red2hdf [options] files                                           */
/* Patchlist:                                                               */
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
#include "hdf.h"
#include "hrc.h"
#include "whipple.h"

/*************************************************************************/
/* command line options                                                  */
/*************************************************************************/
extern int optind;

char *filename;
char *progname;

/*************************************************************************/
/* function declarations                                                 */
/*************************************************************************/
void getoptions();
void usage(void);

int main(int argc, char **argv)
{
  char HDFoutname[FILENAME_MAX];
  FILE *infile;

  struct header head;
  struct frame event;
  double ra, dec;
  double gpstime, livetime, firsttime = 0.0, utcstart;
  int sw_version = 1;
  int nadc = 120;     /* always */
  short zero = 0;

  char *vruninfo_name[] = 
  {"VERSION","NADC","UTC_START","SOURCE_NAME","SOURCE_RA","SOURCE_DEC"};
  char vruninfo_name_str[] = 
    "VERSION,NADC,UTC_START,SOURCE_NAME,SOURCE_RA,SOURCE_DEC";
  char *vevent_name[] = 
  {"CODE","GPSTIME","OSCTIME","LIVETIME","ADC","ELEVATION","AZIMUTH"};
  char vevent_name_str[] = 
    "CODE,GPSTIME,OSCTIME,LIVETIME,ADC,ELEVATION,AZIMUTH";
  int32 file_id, vruninfo_id, vevent_id;
  uint8 *runinfobuf, *eventbuf;

  progname = basename(*argv);
/*************************************************************************/
/* set default command line options                                      */
/*************************************************************************/

/*************************************************************************/
/* get command line options                                              */
/*************************************************************************/
  getoptions(argc, argv);

/****************************************************************************/
/* open reduced file                                                        */
/****************************************************************************/
  filename = *(argv+optind);
  if((infile=fopen(filename,"r")) == NULL){
    fprintf(stderr,"%s: error opening %s\n",progname,*(argv+optind));
    exit(EXIT_FAILURE);
  }

/****************************************************************************/
/* set up HDF file                                                          */
/****************************************************************************/
  sprintf(HDFoutname,"%c%c00%4s.hdf",*filename,*(filename+1),filename+2);
  file_id = Hopen(HDFoutname,DFACC_CREATE,0);
  Vstart(file_id);

  vruninfo_id = VSattach(file_id, -1, "w");
  VSfdefine(vruninfo_id,vruninfo_name[0],DFNT_INT32,   1); /*SW VERSION*/
  VSfdefine(vruninfo_id,vruninfo_name[1],DFNT_INT32,   1); /*NADC      */
  VSfdefine(vruninfo_id,vruninfo_name[2],DFNT_FLOAT64, 1); /*UTC_START */
  VSfdefine(vruninfo_id,vruninfo_name[3],DFNT_CHAR8,  80); /*SOURCE    */
  VSfdefine(vruninfo_id,vruninfo_name[4],DFNT_FLOAT64, 1); /*SOURCE RA */
  VSfdefine(vruninfo_id,vruninfo_name[5],DFNT_FLOAT64, 1); /*SOURCE DEC*/
  VSsetname(vruninfo_id,"10M Run Information");
  VSsetclass(vruninfo_id,"Run Information");
  VSsetfields(vruninfo_id,vruninfo_name_str);

  vevent_id = VSattach(file_id, -1, "w");
  VSfdefine(vevent_id,vevent_name[0],DFNT_INT16,   1); /*EVENT CODE*/
  VSfdefine(vevent_id,vevent_name[1],DFNT_FLOAT64, 1); /*GPSTIME   */
  VSfdefine(vevent_id,vevent_name[2],DFNT_FLOAT64, 1); /*OSCTIME   */
  VSfdefine(vevent_id,vevent_name[3],DFNT_FLOAT64, 1); /*LIVETIME  */
  VSfdefine(vevent_id,vevent_name[4],DFNT_INT16,nadc); /*ADC       */
  VSfdefine(vevent_id,vevent_name[5],DFNT_INT16,   1); /*ELEVATION */
  VSfdefine(vevent_id,vevent_name[6],DFNT_INT16,   1); /*AZIMUTH   */
  VSsetname(vevent_id,"10M Event");
  VSsetclass(vevent_id,"Event");
  VSsetfields(vevent_id,vevent_name_str);

  runinfobuf = (uint8 *)malloc(2*sizeof(int32)+3*sizeof(float64)+
			       80*sizeof(char8));

  printf("%s: reading %s ...\n",progname,filename);
/****************************************************************************/
/* get run information                                                      */
/****************************************************************************/
  head = get_head(infile);
  printf("   runid: %s\n",head.runid);
  printf("   sourename: %s\n",head.sourcename);
  printf("   nevents: %d\n",head.nevents);
  printf("   kdate: %d\n",head.kdate);
  printf("   duration: %lf\n",head.duration);

  memcpy(runinfobuf,&sw_version,4);
  memcpy(runinfobuf+4,&nadc,4);
  utcstart = head.mjd + head.frjd;
  memcpy(runinfobuf+8,&utcstart,8);
  memcpy(runinfobuf+16,head.sourcename,20);
  ra = (double)hhmmss_to_rad((float)head.ra); 
  dec = (double)ddmmss_to_rad((float)head.dec);
  memcpy(runinfobuf+96,&ra,8);
  memcpy(runinfobuf+104,&dec,8);

  eventbuf = (uint8 *)malloc(3*sizeof(float64)+(nadc+3)*sizeof(int16));

  printf("   writing %s ...\n",HDFoutname);
  event = get_frame(infile);
  firsttime =  event.time;
  while(!event.flag){
    if(!event.flag){
      switch(event.code){
      case 8:
	event.code = 6;
	break;
      case 1: case 2:
	event.code = 1;
	break;
      default:
	break;
      }
      memcpy(eventbuf,&event.code,2);
      gpstime = head.mjd + event.time/86400.0L;
      memcpy(eventbuf+2,&gpstime,8);
      memcpy(eventbuf+10,&event.time,8);
      livetime = event.time - firsttime;
      memcpy(eventbuf+18,&livetime,8);
      memcpy(eventbuf+26,event.channel,nadc*2);
      memcpy(eventbuf+nadc*2+26,&zero,2);
      memcpy(eventbuf+nadc*2+28,&zero,2);
      VSwrite(vevent_id, (uchar8 *)eventbuf, 1, FULL_INTERLACE);
    }
    event = get_frame(infile);
  }

  VSwrite(vruninfo_id, (uchar8 *)runinfobuf, 1, FULL_INTERLACE);

  VSdetach(vruninfo_id);
  VSdetach(vevent_id);
  Vend(file_id);
  Hclose(file_id);

  printf("done ...\n");
  return 1;
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
    c = getopt(argc, argv, "h");
    switch(c){
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
  printf("%s [options] filename\n",progname);
  printf(" where [options] are\n");
  printf("   -h     help\n");
  return;
}
