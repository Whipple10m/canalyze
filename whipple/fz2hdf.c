/****************************************************************************/
/* fz2hdf: converts fz files to HDF format                                  */
/* Written by R.W. Lessard 961214                                           */
/* Usage: fz2hdf [options] files                                            */
/* Patchlist:                                                               */
/*            Y2K compliant 12/27/1999                                      */
/*            user definable output 01/28/2000                              */
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
#include "gdf.h"
#include "mbook.h"

#ifdef SUN
#define GDFOPT "SUN"
#endif
#ifdef AXP
#define GDFOPT "AXP"
#endif
#ifdef VMS
#define GDFOPT "VMS"
#endif
#ifdef LINUX
#define GDFOPT "CIO"
#endif

/*************************************************************************/
/* command line options                                                  */
/*************************************************************************/
double osc_freq;
extern int optind;

/*************************************************************************/
/* global variables                                                      */
/*************************************************************************/
char *filename;
char *progname;

/*************************************************************************/
/* diagnostic variables                                                  */
/*************************************************************************/
MHist gpsdt, oscdt, gpsoscdt, rate, rawtime;
MHist trackerr, azimuth, elevation;
char sourceptr[80];
double *hvstatptr;
double *infoptr;
double *faultptr;

/*************************************************************************/
/* function declarations                                                 */
/*************************************************************************/
void getoptions();
void usage(void);
void decode_grs(int *, unsigned *, unsigned *, unsigned *, double, double *, 
		unsigned *, char *);
int nadc;

main(int argc, char **argv)
{
  char HDFoutname[FILENAME_MAX];
  double zero = 0.0;
  int i, ierr = 0, lunit = 10;
  short izero = 0;
  unsigned uval;

  char *vruninfo_name[] = 
  {"VERSION","NADC","UTC_START","SOURCE_NAME","SOURCE_RA","SOURCE_DEC"};
  char vruninfo_name_str[] = 
    "VERSION,NADC,UTC_START,SOURCE_NAME,SOURCE_RA,SOURCE_DEC";
  char *vevent_name[] = 
  {"CODE","GPSTIME","OSCTIME","LIVETIME","ADC","ELEVATION","AZIMUTH"};
  char vevent_name_str[] = 
    "CODE,GPSTIME,OSCTIME,LIVETIME,ADC,ELEVATION,AZIMUTH";
  int16 *code, *hval;
  int32 file_id, vruninfo_id, vevent_id;
  uint8 *runinfobuf, *eventbuf;

  double utcstart, utcend;
  double gpsfrday, firstgps, secondgps, lastgps;
  double osctime, avoscfreq = 0.0, lastosc;
  double livetime, add_live_time;
  int version, /*nadc,*/  nminutes;
  int gpssecond, firstsecond, secondsecond, lastsecond;
  int scalersecond, after_first_gpssecond;
  int prev_live_sec, prev_live_ns;

  int missedmarkers = 0;
  unsigned long long scaler, firstscaler, lastscaler = 0;
  unsigned mark, firstmark, secondmark, lastmark, lastmissedmark;
  unsigned gateclose, lastgateclose;
  double twopow32 = (double)~0U + 1.0L;
  unsigned wrappedgateclose, wrappedscaler, wrappedmark;
  
  char qual;
  double secs;
  unsigned dayofyear, hours, mins, istatus;

  mbook1(&gpsdt,500,0.0,500.0);        /* ms  */
  mbook1(&oscdt,500,0.0,500.0);        /* ms  */
  mbook1(&gpsoscdt,1000,-500.0,500.0); /* ms  */
  mbook1(&rawtime,524288,0.0,1680.0);  /* s   */
  mbook1(&trackerr,50,0.0,1.0);        /* deg */
  mbook1(&azimuth,360,-180.0,180.0);   /* deg */
  mbook1(&elevation,90,0.0,90.0);      /* deg */

  progname = basename(*argv);
/*************************************************************************/
/* set default command line options                                      */
/*************************************************************************/
  osc_freq = 0.0;

/*************************************************************************/
/* get command line options                                              */
/*************************************************************************/
  getoptions(argc, argv);

  uval = 1; gdf_option_("CIO",&uval,&ierr,strlen("CIO"));
  uval = 1; gdf_option_("RESET",&uval,&ierr,strlen("RESET"));

  gdf_init_("Z",&ierr,1);
  if(ierr){
    fprintf(stderr,"%s: error initializing GDF\n",progname);
    exit(EXIT_FAILURE);
  }
  filename = *(argv+optind);

/****************************************************************************/
/* get run information                                                      */
/****************************************************************************/
  gdf_open_(&lunit,filename,"R",&ierr,strlen(filename),1);
  if(ierr){
    fprintf(stderr,"%s: error opening %s\n",progname,filename);
    exit(EXIT_FAILURE);
  }
  printf("%s: reading %s ...\n",progname,filename);
  printf("   searching for run info ...\n");
  version = 0; nadc = 0;
  while((ierr != 1) && ((version == 0) || (nadc == 0))){
    gdf_read_(&lunit," ",&ierr,1);
    if(gdf_data_.gdf_run.new && (ierr==0)){
      version = gdf_data_.gdf_run.version;
      utcstart = gdf_data_.gdf_run.utc_start;
      utcend = gdf_data_.gdf_run.utc_end;
      nminutes = (int)((utcend - (int)utcend)*1440.0) - 
		       (int)((utcstart - (int)utcstart)*1440.0);
      printf("      run duration: %d mins\n",nminutes);
    }
    if(gdf_data_.gdf_ev10.new && (ierr == 0) && (nadc == 0)){
      nadc = gdf_data_.gdf_ev10.nadc;
      printf("      no. channels: %d\n",nadc);
    }
  }
  printf("   done ...\n");
  gdf_close_(&lunit,&ierr);
  if(version == 0) printf("   **warning** no run info\n");
  if(nadc == 0){
    printf("   **error** no events\ndone ...\n");
    exit(EXIT_FAILURE);
  }

  mbook1(&rate,nminutes,0.0,(double)nminutes);

  hvstatptr = (double *) malloc(nadc*sizeof(double));
  for(i=0;i < nadc;i++){
    *(hvstatptr+i) = 1.0;
  }

  faultptr = (double *) malloc(6*sizeof(double));

/****************************************************************************/
/* oscillator calibration                                                   */
/****************************************************************************/
  gdf_open_(&lunit,filename,"R",&ierr,strlen(filename),1);
  if(osc_freq == 0.0){
    printf("   calibrating oscillator ...\n");
    *(faultptr) = 1.0;
  }else{
    printf("   looking for GPS second markers ...\n");
    *(faultptr) = 0.0;
  }
  firstmark = 0; firstscaler = 0ULL; wrappedmark = 0; wrappedscaler = 0;
  firstgps = 0.0;
  while(ierr != 1){
    gdf_read_(&lunit," ",&ierr,1);
    if(gdf_data_.gdf_ev10.new && (ierr == 0)){
/****************************************************************************/
/* oscillator calibration - pre version 76                                  */
/****************************************************************************/
      gpsfrday = gdf_data_.gdf_ev10.utc-
	(double)((int)gdf_data_.gdf_ev10.utc);
      gpssecond = (int)(gpsfrday*86400.0);
      if(version < 76){
	mark = gdf_data_.gdf_ev10.mark_gps;
	if(firstmark == 0){
	  firstmark = mark; lastmark = firstmark;
	}
/****************************************************************************/
/* new GPS second marker                                                    */
/****************************************************************************/
	if(mark != lastmark){
	  if(lastmark == firstmark){
/****************************************************************************/
/* first mark                                                               */
/****************************************************************************/
	    secondmark = mark;
	    secondgps = (double)gpssecond;
	    secondsecond = gpssecond;
	    printf("      SECOND MARK %u SECOND GPS %lf [s]\n",
		   secondmark,secondgps);	  
	    if(osc_freq != 0.0)
	      break;
	  }
	  if(mark < lastmark)
	    ++wrappedmark;
	  lastmark = mark;
	  lastgps = (double)gpssecond;
	}
      }else{
/****************************************************************************/
/* oscillator calibration - post version 76                                 */
/****************************************************************************/
	scaler = 
	  (unsigned long long)(gdf_data_.gdf_ev10.elapsed_sec)*1000000000 + 
	  (unsigned long long)(gdf_data_.gdf_ev10.elapsed_ns);
/****************************************************************************/
/* new GPS second marker                                                    */
/****************************************************************************/
	if(gdf_data_.gdf_ev10.trigger & 1){
	  if((firstscaler == 0) && (firstgps == 0.0)){
	    firstscaler = scaler;
	    firstgps = (double)gpssecond;
	    firstsecond = gpssecond;
	    printf("      FIRST SCALER %llu FIRST GPS %lf [s]\n",
		   firstscaler,firstgps);
	    if(osc_freq != 0.0)
	      break;
	  }
	  if(scaler < lastscaler){
	    ++wrappedscaler;
	  }
	  lastscaler = scaler;
	  lastgps = (double)gpssecond;
	}
      }
    }
  }
  
  if(version < 76){
    if(osc_freq == 0.0){
      printf("      LAST   MARK %u WRAPS %u LAST GPS %lf [s]\n",
	     lastmark,wrappedmark,lastgps);
      avoscfreq = ((double)lastmark + (double)wrappedmark*twopow32 - 
		   (double)secondmark)/(lastgps-secondgps);
      printf("      AVERAGE OSCILLATOR FREQUENCY %lf [Hz]\n",avoscfreq);
    }else{
      avoscfreq = osc_freq*1000000.0;
    }
  }else{
    if(osc_freq == 0.0){
      printf("      LAST  SCALER %llu WRAPS %u LAST GPS %lf [s]\n",
	     lastscaler,wrappedscaler,lastgps);
      if(lastgps == firstgps){
	avoscfreq = 1.0e9;
	printf("      **warning** calibration impossible, using 10^9[ns/s]\n");
	*(faultptr) = 0.0;
      }else{
	avoscfreq = ((double)lastscaler + (double)wrappedscaler*twopow32 - 
		     (double)firstscaler)/(lastgps-firstgps);
	printf("      AVERAGE OSCILLATOR FREQUENCY %lf [ns/s]\n",avoscfreq);
      }
    }else{
      avoscfreq = 1.0e9;
    }
  }

  *(faultptr+1) = avoscfreq;
  gdf_close_(&lunit,&ierr);
  printf("   done ...\n");

  gdf_open_(&lunit,filename,"R",&ierr,strlen(filename),1);
/****************************************************************************/
/* set up HDF file                                                          */
/****************************************************************************/
  strncpy(HDFoutname,(char *)basename(filename),8); HDFoutname[8]='\0';
  strcat(HDFoutname,".hdf");

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
  memcpy(runinfobuf,&version,4);
  memcpy(runinfobuf+4,&nadc,4);
  memcpy(runinfobuf+8,&utcstart,8);
  memcpy(runinfobuf+16,"UNKNOWN",80);
  memcpy(runinfobuf+96,&zero,8);
  memcpy(runinfobuf+104,&zero,8);

  eventbuf = (uint8 *)malloc(3*sizeof(float64)+(nadc+3)*sizeof(int16));
  code = (int16 *)malloc(sizeof(int16));
  hval = (int16 *)malloc(sizeof(int16)*nadc);

  printf("   writing %s ...\n",HDFoutname);
  prev_live_sec = 0; prev_live_ns = 0; add_live_time = 0.0;
  after_first_gpssecond = 0; 
  firstgps = 0.0; lastgps = 0.0; lastosc = 0.0;
  lastmark = 0; lastmissedmark = 0; lastgateclose = 0; lastsecond = 0;
  wrappedgateclose = 0; wrappedmark = 0; 
  *(faultptr+2) = 0.0; *(faultptr+3) = 0.0; 
  *(faultptr+4) = 0.0; *(faultptr+5) = 0.0;
  while(ierr != 1){
    gdf_read_(&lunit," ",&ierr,1);
    if(gdf_data_.gdf_track[0].new && (gdf_data_.gdf_track[0].utc > utcstart) 
       && (ierr == 0)){
      memcpy(runinfobuf+16,gdf_data_.gdf_track[0].source,80);
      memcpy(runinfobuf+96,&(gdf_data_.gdf_track[0].rasc_today),8);
      memcpy(runinfobuf+104,&(gdf_data_.gdf_track[0].decl_today),8);
      mfill(&trackerr,gdf_data_.gdf_track[0].deviation*180/M_PI,
	    gdf_data_.gdf_track[0].deviation*180/M_PI);
      mfill(&azimuth,gdf_data_.gdf_track[0].azimuth*180/M_PI,
	    gdf_data_.gdf_track[0].azimuth*180/M_PI);
      mfill(&elevation,gdf_data_.gdf_track[0].elevation*180/M_PI,
	    gdf_data_.gdf_track[0].elevation*180/M_PI);

    }
    if(gdf_data_.gdf_hv[0].new && (gdf_data_.gdf_hv[0].utc > utcstart) && 
       (ierr == 0)){
      for(i=0;i < nadc;i++)
	if(!(gdf_data_.gdf_hv[0].status[i] & 1))
	  *(hvstatptr+i) = 0.0;
    }
    if(gdf_data_.gdf_ev10.new && (ierr == 0)){
      gpsfrday = gdf_data_.gdf_ev10.utc-(double)((int)gdf_data_.gdf_ev10.utc);
      gpssecond = (int)(gpsfrday*86400.0L);
      if(firstgps == 0.0)
	firstgps = gpsfrday;
/****************************************************************************/
/* oscillator times - pre version 76                                        */
/****************************************************************************/
      if(version < 76){
/****************************************************************************/
/* scaler at last GPS second marker                                         */
/****************************************************************************/
	mark = gdf_data_.gdf_ev10.mark_gps;
	gateclose = gdf_data_.gdf_ev10.gate_close;	
	/*	printf("MARK %u LMARK %u GATE %u GPS %d LGPS %d\n",
	  mark,lastmark,gateclose,gpssecond,lastsecond); */
/****************************************************************************/
/* count scaler wraps                                                       */
/****************************************************************************/
	if(gateclose < lastgateclose) ++wrappedgateclose;
	if(mark < lastmark) ++wrappedmark;
/****************************************************************************/
/* first marker is bogus so use second marker to calc osc times             */
/****************************************************************************/
	if(mark == firstmark){
	  osctime = (double)gpssecond + (double)(secondsecond - gpssecond) -
	    ((double)(secondmark-gateclose))/avoscfreq;
	}else{
/****************************************************************************/
/* GPS second marker is latched on new GPS second after gateclose           */
/****************************************************************************/
	  if((gateclose < mark) && (wrappedgateclose == wrappedmark)){
	    printf("      **info** second marker arrived after gateclose\n"); 
	    mark = lastmark;
	  }
/****************************************************************************/
/* old mark - new gps second (missed marker) -> find matching GPS second    */
/****************************************************************************/
	  missedmarkers = gpssecond - lastsecond -
	    (int)(((double)(mark - lastmark)/avoscfreq) + 0.5);
	  if((missedmarkers > 0) && (lastmark != firstmark)){
	    if(mark != lastmissedmark)
	      printf("      **info** missed %d second markers\n",
		     missedmarkers);
	    lastmissedmark = mark;
	    gpssecond -= missedmarkers;
	    *(faultptr+2) += missedmarkers;
	  }
/****************************************************************************/
/* new mark - old gps second -> use last mark                               */
/****************************************************************************/
	  if((mark != lastmark) && (gpssecond == lastsecond)){
	    printf("      **info** marker arrived before GPS second\n");
	    mark = lastmark;
	  }
	  osctime = (double)gpssecond + (double)(gateclose - mark)/avoscfreq;
	}
	lastsecond = gpssecond;	lastmark = mark; lastgateclose = gateclose;
      }else{
/****************************************************************************/
/* oscillator times - post version 76                                       */
/****************************************************************************/
	scaler = 
	  (unsigned long long)(gdf_data_.gdf_ev10.elapsed_sec)*
	  1000000000ULL + 
	  (unsigned long long)(gdf_data_.gdf_ev10.elapsed_ns);
	if(gdf_data_.gdf_ev10.trigger & 1){
	  lastscaler = scaler;
	  if(gpssecond == lastsecond)            /* temporary DACQ bug fix */
	    scalersecond = gpssecond + 1;
	  else if(gpssecond = (lastsecond - 1 )) /* temporary DACQ bug fix */
	    scalersecond = gpssecond + 2;
	  else
	    scalersecond = gpssecond;
	  after_first_gpssecond = 1;
	}
	if(after_first_gpssecond)
	  osctime = (double)scalersecond + 
	    (double)((scaler - lastscaler))/avoscfreq;
	else
	  osctime = (double)gpssecond + (double)(firstsecond - gpssecond) -
	    ((double)(firstscaler - scaler))/avoscfreq;
	lastsecond = gpssecond;
      }
/****************************************************************************/
/* live time                                                                */
/****************************************************************************/
      if(gdf_data_.gdf_ev10.live_sec < prev_live_sec){
	if(gdf_data_.gdf_ev10.live_sec == 0){
	  printf("      **info** detected live time reset at %.3lfs\n",
		 livetime);
	  printf("               this counter %ds %dns, last %ds %dns\n",
		 gdf_data_.gdf_ev10.live_sec,gdf_data_.gdf_ev10.live_ns,
		 prev_live_sec,prev_live_ns);
	}else{ 
	  printf("      **error** found live time counter fault at %.3lfs\n",
		 livetime);
	  printf("                assuming a reset + time gap\n");
	  printf("                this counter %ds %dns",
		 gdf_data_.gdf_ev10.live_sec,
		 gdf_data_.gdf_ev10.live_ns);
	  printf(" last %ds %dns\n",prev_live_sec,prev_live_ns);
	}
	add_live_time += ((double)prev_live_sec + (double)prev_live_ns/1.0e9);
      }
      prev_live_sec = gdf_data_.gdf_ev10.live_sec; 
      prev_live_ns = gdf_data_.gdf_ev10.live_ns;
      livetime = (double)gdf_data_.gdf_ev10.live_sec + 
	(double)gdf_data_.gdf_ev10.live_ns/1.0e9 + add_live_time;
#ifdef DEBUG
      printf("SCALER %llu LASTSCALER %llu SCALERSEC %d\n",scaler,lastscaler,
	     scalersecond);
      printf("UTC %lf ELAPSED %d %d LIVE %d %d\n",gdf_data_.gdf_ev10.utc,
	     gdf_data_.gdf_ev10.elapsed_sec,gdf_data_.gdf_ev10.elapsed_ns,
	     gdf_data_.gdf_ev10.live_sec,gdf_data_.gdf_ev10.live_ns); 
      printf("GPS %lf OSC %lf EL %lf LV %lf\n",gpsfrday*86400.0, 
	     osctime,(double)gdf_data_.gdf_ev10.elapsed_sec + 
	     (double)gdf_data_.gdf_ev10.elapsed_ns/1.0e9,livetime);
#endif
/****************************************************************************/
/* check GPS status bits                                                    */
/****************************************************************************/
      if(version > 76){
	decode_grs(gdf_data_.gdf_ev10.grs_clock,&dayofyear,&hours,&mins,
		   avoscfreq,&secs,&istatus,&qual);
/****************************************************************************/
/* okay done ... now for some moaning                                       */
/****************************************************************************/
	if(!(istatus & 1)){
	  printf("      **warning** GRS2: no 1 pps\n");
	  *(faultptr+2) = 1.0;
	}
	if(!(istatus & 2)){
	  printf("      **warning** GRS2: 10MHz clock unsteady\n");
	  *(faultptr+2) = 1.0;
	}
	if(!(istatus & 4)){
	  printf("      **warning** GRS2: RS-232 data invalid and not incrementing\n");
	  *(faultptr+2) = 1.0;
	}
	if(!(istatus & 8)){
	  printf("      **warning** GRS2: RS-232 data does not match internal time\n");
	  *(faultptr+2) = 1.0;
	}
	/*if(istatus == 15){
	  printf("GRS2 returned: %u days, %2u:%2u:%9.6lf hh.mm.ss TQC %c\n",
	  dayofyear,hours,mins,secs,qual); */
      }
      if(osctime < lastosc){
	printf("      **error** osc time: %lf [s] previous %lf [s]\n",
	       osctime,lastosc);
	*(faultptr+4) += 1.0;
      }
      if(gpsfrday < lastgps){
	printf("      **error** gps time: %lf [s] previous %lf [s]\n",
	       gpsfrday*86400.0,lastgps*86400.0);
	*(faultptr+3) += 1.0;
      }
      if(fabs((gpsfrday*86400.0 - osctime)*1000.0) > 1.0){
	printf("      **warning** osc/gps time: %.4lf/%.4lf diff %.2lf [ms]\n",
	       osctime,gpsfrday*86400.0,
	       fabs((gpsfrday*86400.0 - osctime)*1000.0));
	*(faultptr+5) += 1.0;
      }

      if(lastgps != 0.0){
	mfill(&gpsdt,(gpsfrday - lastgps)*86400.0*1000.0,
	      (gpsfrday - lastgps)*86400.0*1000.0);
	mfill(&oscdt,(osctime - lastosc)*1000.0,(osctime - lastosc)*1000.0);
      }
      mfill(&gpsoscdt,(gpsfrday*86400.0 - osctime)*1000.0,
	    (gpsfrday*86400.0 - osctime)*1000.0);
      lastgps = gpsfrday;
      lastosc = osctime;

      if(version < 76)
	*code = 4;
      else
	*code = (int16)gdf_data_.gdf_ev10.trigger;
      if(*code & 4){
	mfill(&rawtime,(gpsfrday - firstgps)*86400.0,
	      (gpsfrday - firstgps)*86400.0);
	mfill(&rate,livetime/60.0,livetime/60.0);
      }
      memcpy(eventbuf,code,2);
      memcpy(eventbuf+2,&(gdf_data_.gdf_ev10.utc),8);
      memcpy(eventbuf+10,&osctime,8);
      memcpy(eventbuf+18,&livetime,8);
#ifdef SUN
      for(i=0;i < nadc;i++)
	*(hval+i) = (int16)(gdf_data_.gdf_ev10.adc[i]);
      memcpy(eventbuf+26,hval,nadc*2);
      *hval = (int16)(gdf_data_.gdf_ev10.track[0]);
      memcpy(eventbuf+nadc*2+26,hval,2);
      *hval = (int16)(gdf_data_.gdf_ev10.track[1]);
      memcpy(eventbuf+nadc*2+28,hval,2);
#else
      memcpy(eventbuf+26,gdf_data_.gdf_ev10.adc,nadc*2);
      memcpy(eventbuf+nadc*2+26,&(gdf_data_.gdf_ev10.track[0]),2);
      memcpy(eventbuf+nadc*2+28,&(gdf_data_.gdf_ev10.track[1]),2);
#endif
      VSwrite(vevent_id, (uchar8 *)eventbuf, 1, FULL_INTERLACE);
    }
    if(gdf_data_.gdf_fr10.new && (ierr == 0) && (version < 76)){
      *code = 1;
      memcpy(eventbuf,code,2);
      memcpy(eventbuf+2,&(gdf_data_.gdf_fr10.utc),8);
      memcpy(eventbuf+10,&osctime,8);
      memcpy(eventbuf+18,&livetime,8);
#ifdef SUN
      for(i=0;i < nadc;i++)
	*(hval+i) = (int16)(gdf_data_.gdf_fr10.ped_adc1[i]);
      memcpy(eventbuf+26,hval,nadc*2);
#else
      memcpy(eventbuf+26,gdf_data_.gdf_fr10.ped_adc1,nadc*2);
#endif
      memcpy(eventbuf+nadc*2+26,&izero,2);
      memcpy(eventbuf+nadc*2+28,&izero,2);
      VSwrite(vevent_id, (uchar8 *)eventbuf, 1, FULL_INTERLACE);

#ifdef SUN
      for(i=0;i < nadc;i++)
	*(hval+i) = (int16)(gdf_data_.gdf_fr10.ped_adc2[i]);
      memcpy(eventbuf+26,hval,nadc*2);
#else
      memcpy(eventbuf+26,gdf_data_.gdf_fr10.ped_adc2,nadc*2);
#endif
      VSwrite(vevent_id, (uchar8 *)eventbuf, 1, FULL_INTERLACE);
    }
  }
  VSwrite(vruninfo_id, (uchar8 *)runinfobuf, 1, FULL_INTERLACE);

  VSdetach(vruninfo_id);
  VSdetach(vevent_id);
  Vend(file_id);
  Hclose(file_id);
  gdf_close_(&lunit,&ierr);
  printf("   done ...\n");
  printf("   run live time %lf mins\n",livetime/60.0);

  infoptr = (double *) malloc(8*sizeof(double));

  *(infoptr) = (double)nminutes;
  *(infoptr+1) = (double)nadc;
  *(infoptr+2) = utcstart;
  *(infoptr+3) = livetime/60.0;
  *(infoptr+4) = gdf_data_.gdf_track[0].rasc_today; 
  *(infoptr+5) = gdf_data_.gdf_track[0].decl_today; 
  memcpy(sourceptr,runinfobuf+16,80);

#ifdef HAVE_MATLAB
  MATLAB_output();
#endif
#ifdef HAVE_USER_FZ2HDF_OUTPUT
  USER_FZ2HDF_output();
#endif
  printf("done ...\n");

  gdf_exit_(&ierr);

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
    c = getopt(argc, argv, "ho:");
    switch(c){
    case 'h':
      usage();
      exit(EXIT_SUCCESS);
    case 'o':
      if(sscanf(optarg,"%lf",&osc_freq) != 1)
	printf("%s: invalid option\n",progname);
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
  printf("%s [options] filename.fz\n",progname);
  printf(" where [options] are\n");
  printf("   -h     help\n");
  printf("   -o %%f  oscillator freq [MHz]- no calibration\n");
  return;
}

/*************************************************************************/
/* function: decode_grs                                                  */
/* passed:   grstime bits                                                */
/* returns:  decoded time                                                */
/* purpose:  decode GRS GPS clock interface                              */
/*************************************************************************/
void decode_grs(int *grstime, unsigned *days, 
		unsigned *hours, unsigned *mins, 
		double oscfreq, double *secs,
		unsigned *status, char *qual)
{
  unsigned int seconds, clock_tics;

  clock_tics = grstime[0];
  seconds = (grstime[1] & 15) + 10*((grstime[1] >> 4) & 15);
  *secs = (double)seconds + (double)clock_tics/oscfreq;

  *mins = ((grstime[1] >> 8) & 15) + 10*((grstime[1] >> 12) & 15);
  *hours = ((grstime[1] >> 16) & 15) + 10*((grstime[1] >> 20) & 15);

  *days = (grstime[2] & 15) +
    10*((grstime[2] >> 4) & 15) +
      100*((grstime[2] >> 8) & 15);
  
  *qual = (char)(((grstime[2] >> 12) & 15));
  *status = ((grstime[2] >> 16) & 15);

  return;
}

#ifdef HAVE_MATLAB
#include "mat.h"

int MATLAB_output(void)
{
  char MAToutname[FILENAME_MAX];
  mxArray *mhvstat;
  mxArray *msource, *minfo, *mfault;
  mxArray *mgpsdt, *moscdt, *mgpsoscdt, *mrawtime, *mrate;
  mxArray *mtrackerr, *mazimuth, *melevation;

  MATFile *matfile;

  strncpy(MAToutname,(char *)basename(filename),8); MAToutname[8]='\0';
  strcat(MAToutname,"d.mat");
  printf("   writing MATLAB diagnostics file %s ...\n",MAToutname);

  // SJF - this next macro copies the memory beloning to a MATLAB matrix
  // into another MATLAB matrix. This must be done since MATLAB wants to
  // free both of them when they are finished -- you cannot just assign
  // the pointer as RWL had been doing!

  // 1-d matrix (vector)  COPYMBOOK3(mtime,output_results.photon.time,100000);

#define COPYMBOOK1(A,M) memcpy(mxGetPr(A),(M).matrixptr,(M).nibins*sizeof(double))
#define COPYMBOOK3(A,M,N) memcpy(mxGetPr(A),M,N*sizeof(double))

  mhvstat = mxCreateDoubleMatrix((int)(*(infoptr+1)),1,mxREAL);
  mxSetName(mhvstat,"hvstat"); // mxSetPr(mhvstat,hvstatptr);
  COPYMBOOK3(mhvstat,hvstatptr,nadc);

  msource = mxCreateString(sourceptr);
  mxSetName(msource,"source");

  minfo = mxCreateDoubleMatrix(8,1,mxREAL);
  mxSetName(minfo,"minfo"); // mxSetPr(minfo,infoptr);
  COPYMBOOK3(minfo,infoptr,8);

  mfault = mxCreateDoubleMatrix(6,1,mxREAL);
  mxSetName(mfault,"mfault"); // mxSetPr(mfault,faultptr);
  COPYMBOOK3(mfault,faultptr,6);


  mgpsdt = mxCreateDoubleMatrix(gpsdt.nibins,1,mxREAL);
  mxSetName(mgpsdt,"gpsdt"); // SJF WRONG: mxSetPr(mgpsdt,gpsdt.matrixptr);
  COPYMBOOK1(mgpsdt,gpsdt);

  moscdt = mxCreateDoubleMatrix(oscdt.nibins,1,mxREAL);
  mxSetName(moscdt,"oscdt"); // SJF WRONG: mxSetPr(moscdt,oscdt.matrixptr);
  COPYMBOOK1(moscdt,oscdt);

  mgpsoscdt = mxCreateDoubleMatrix(gpsoscdt.nibins,1,mxREAL);
  mxSetName(mgpsoscdt,"gpsoscdt"); // SJF WRONG: mxSetPr(mgpsoscdt,gpsoscdt.matrixptr);
  COPYMBOOK1(mgpsoscdt,gpsoscdt);

  mrawtime = mxCreateDoubleMatrix(rawtime.nibins,1,mxREAL);
  mxSetName(mrawtime,"rawtime"); // SJF WRONG: mxSetPr(mrawtime,rawtime.matrixptr);
  COPYMBOOK1(mrawtime,rawtime); 

  mrate = mxCreateDoubleMatrix(rate.nibins,1,mxREAL);
  mxSetName(mrate,"rate"); // SJF WRONG: mxSetPr(mrate,rate.matrixptr);
  COPYMBOOK1(mrate,rate);

  mtrackerr = mxCreateDoubleMatrix(trackerr.nibins,1,mxREAL);
  mxSetName(mtrackerr,"trackerr"); // SJF WRONG: mxSetPr(mtrackerr,trackerr.matrixptr);
  COPYMBOOK1(mtrackerr,trackerr);

  mazimuth = mxCreateDoubleMatrix(azimuth.nibins,1,mxREAL);
  mxSetName(mazimuth,"azi"); // SJF WRONG: mxSetPr(mazimuth,azimuth.matrixptr);
  // JPF I had to change the matlab variable name from "azimuth" to "azi"
  COPYMBOOK1(mazimuth,azimuth);

  melevation = mxCreateDoubleMatrix(elevation.nibins,1,mxREAL);
  mxSetName(melevation,"elev"); // SJF WRONG: mxSetPr(melevation,elevation.matrixptr);
  // JPF I had to change the matlab variable name from "elevation" to "elev"
  COPYMBOOK1(melevation,elevation);

  /*  mgpsdt = mxCreateDoubleMatrix(gpsdt.nibins,1,mxREAL);
      mxSetName(mgpsdt,"gpsdt"); mxSetPr(mgpsdt,gpsdt.matrixptr);

      moscdt = mxCreateDoubleMatrix(oscdt.nibins,1,mxREAL);
      mxSetName(moscdt,"oscdt"); mxSetPr(moscdt,oscdt.matrixptr);

      mgpsoscdt = mxCreateDoubleMatrix(gpsoscdt.nibins,1,mxREAL);
      mxSetName(mgpsoscdt,"gpsoscdt"); mxSetPr(mgpsoscdt,gpsoscdt.matrixptr);

      mrawtime = mxCreateDoubleMatrix(rawtime.nibins,1,mxREAL);
      mxSetName(mrawtime,"rawtime"); mxSetPr(mrawtime,rawtime.matrixptr);

      mrate = mxCreateDoubleMatrix(rate.nibins,1,mxREAL);
      mxSetName(mrate,"rate"); mxSetPr(mrate,rate.matrixptr);

      mtrackerr = mxCreateDoubleMatrix(trackerr.nibins,1,mxREAL);
      mxSetName(mtrackerr,"trackerr"); mxSetPr(mtrackerr,trackerr.matrixptr);

      mazimuth = mxCreateDoubleMatrix(azimuth.nibins,1,mxREAL);
      mxSetName(mazimuth,"azimuth"); mxSetPr(mazimuth,azimuth.matrixptr);

      melevation = mxCreateDoubleMatrix(elevation.nibins,1,mxREAL);
      mxSetName(melevation,"elevation"); mxSetPr(melevation,elevation.matrixptr);*/

  matfile = matOpen(MAToutname,"w");
  matPutArray(matfile,msource); matPutArray(matfile,minfo);
  matPutArray(matfile,mgpsdt); matPutArray(matfile,moscdt);
  matPutArray(matfile,mgpsoscdt); matPutArray(matfile,mrawtime);
  matPutArray(matfile,mtrackerr); matPutArray(matfile,mazimuth);
  matPutArray(matfile,melevation); matPutArray(matfile,mrate);
  matPutArray(matfile,mhvstat); matPutArray(matfile,mfault);
  matClose(matfile);

  mxDestroyArray(msource); mxDestroyArray(minfo);
  mxDestroyArray(mgpsdt); mxDestroyArray(moscdt);
  mxDestroyArray(mgpsoscdt); mxDestroyArray(mrawtime);
  mxDestroyArray(mtrackerr); mxDestroyArray(mazimuth);
  mxDestroyArray(melevation); mxDestroyArray(mrate);
  mxDestroyArray(mhvstat);

  return 1;
}
#endif





