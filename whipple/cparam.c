/***************************************************************************/
/* cparam: parameterizes hdf files                                         */
/* Written by R.W. Lessard 960606                                          */
/* Usage: cparam [options] *.hdf                                           */
/* Patchlist:                                                              */
/* Version: HDF                                                            */
/*          Y2K compliant 12/27/1999                                       */
/*          MPI 02/24/2000
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
#include "whipple.h"
#include "sla.h"
#include "nr.h"
#ifdef MPI
#include "mpi.h"
#endif

/*************************************************************************/
/* defines                                                               */
/*************************************************************************/
#define DEFAULT_PICTHRESH 4.25
#define DEFAULT_BNDTHRESH 2.25
#define DEFAULT_EFACTOR   1.37
#define DEFAULT_PATH      "./"
#define SINLATITUDE       0.5252241379
#define COSLATITUDE       0.8509639269
#define WHIPLONG          1.935190
#define WHIPLAT           0.552978

/*************************************************************************/
/* global variables                                                      */
/*************************************************************************/
char *filename;
char *progname;
int myid = 0, master = 0;

/*************************************************************************/
/* position of pmts                                                      */
/*************************************************************************/
struct coords pmt;
int **pmt_neighbours;
float pmt_spacing;

/*************************************************************************/
/* param information                                                     */
/*************************************************************************/
struct parinfo_s{
  char source_name[80];
  char gainlabel[3];
  char n2id[9];
  char padlabel[4];
  char pairid[9];
  double utc_start;
  double source_ra;
  double source_dec;
  float picture_threshold;
  float boundary_threshold;
  float efactor;
  float xo;
  float yo;
  int version;
  int nadc;
  int ntrig;
  int npmt;
  int idate;
  int itime;
  int  ntubesoff;
  int  *tubesoff;
} parinfo;

/*************************************************************************/
/* hillas parameters                                                     */
/*************************************************************************/
struct par_s {
  double gpstime;
  double osctime;
  double livetime;
  float length;
  float width;
  float miss;
  float distance;
  float azwidth;
  float frac[3];
  float maximum[3];
  float asymm;
  float xc;
  float yc;
  float xo;
  float yo;
  int size;
  int location[3];
};

/*************************************************************************/
/* command line options                                                  */
/*************************************************************************/
char path_to_db[FILENAME_MAX];
int fixed;
int pad;
int pairid;
int n2id;
int verbose;
int wavelet;
extern int optind;

/*************************************************************************/
/* function declarations                                                 */
/*************************************************************************/
int cimpad();
int cimclean();
int cimgain();
int get_hillas_parameters(float *, struct par_s *, float, float, int);
void getmax();
void getfrac(float *, float *, float);
void getoptions(int argc, char **argv);
     /*void getoptions();*/
void usage();
void derot(float *, float *, double);
double derot_angle(double);

int main(int argc, char **argv)
{
  char fileid[9];
  char HDF_outname[FILENAME_MAX];
  float *peds, *pairpeds, *pedvars, *pairpedvars, *maxpedvars;
  float *gains, *fadc;
  char cr, cd;
  double ra2000, dec2000, j2000;
  double dfr, epoch, siderealtime;
  int iyy, imm, idd, nyy, ndd;
  int i, j, k, l, itime[4], ra[4], dec[4];
  int *pairtubesoff;
  int ntubesoff, *picture, *boundary;
  short elevation, azimuth;
  short code, *adc;

  struct par_s par;

  int32 infile_id, vevent_ref, vevent_id, vruninfo_ref, vruninfo_id;
  uchar8 *runinfobuf, *eventbuf;

  int32 outfile_id, vparinfo_ref, vparinfo_id, vpar_ref, vpar_id;
  char *vparinfo_name[] = {"VERSION","NADC","UTC_START","IDATE","ITIME",
			   "SOURCE_NAME","SOURCE_RA","SOURCE_DEC",
			   "GAINLABEL","N2ID","PADLABEL","NTUBESOFF",
			   "TUBESOFF","PICTURE","BOUNDARY","EFACTOR"};
  char vparinfo_name_str[] = "VERSION,NADC,UTC_START,IDATE,ITIME,SOURCE_NAME,SOURCE_RA,SOURCE_DEC,GAINLABEL,N2ID,PADLABEL,NTUBESOFF,TUBESOFF,PICTURE,BOUNDARY,EFACTOR";
  char *vpar_name[] = {"GPSTIME","OSCTIME","LIVETIME","LENGTH","WIDTH","MISS",
			 "DISTANCE","AZWIDTH","FRAC","SIZE","LOC","MAX",
			 "ASYMM","XC","YC","XO","YO"};
  char vpar_name_str[] = "GPSTIME,OSCTIME,LIVETIME,LENGTH,WIDTH,MISS,DISTANCE,AZWIDTH,FRAC,SIZE,LOC,MAX,ASYMM,XC,YC,XO,YO";
  uint8 *parinfobuf, **parbuf;

  int numproc = 1, ierr = 0, proc = 1;

#ifdef MPI
/*************************************************************************/
/* initialize MPI                                                        */
/*************************************************************************/
  MPI_Status *status;
  MPI_Request *requests;

  int parinfo_blockcounts[4] = {105,3,5,8};
  MPI_Datatype parinfo_types[4] = {MPI_CHAR,MPI_DOUBLE,MPI_FLOAT,MPI_INT};
  MPI_Aint parinfo_displs[4];
  MPI_Datatype parinfo_datatype;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numproc);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);

/*************************************************************************/
/* define some MPI data types                                            */
/*************************************************************************/

  MPI_Address(&parinfo.source_name,&parinfo_displs[0]);
  MPI_Address(&parinfo.utc_start,&parinfo_displs[1]);
  MPI_Address(&parinfo.picture_threshold,&parinfo_displs[2]);
  MPI_Address(&parinfo.version,&parinfo_displs[3]);
  for(i=3;i >= 0; i--)
    parinfo_displs[i] -= parinfo_displs[0];
  MPI_Type_struct(4,parinfo_blockcounts,parinfo_displs,
		  parinfo_types,&parinfo_datatype);
  MPI_Type_commit(&parinfo_datatype);

  if(numproc > 1){
    requests = (MPI_Request *)malloc(sizeof(MPI_Request)*(numproc-1));
    status = (MPI_Status *)malloc(sizeof(MPI_Status)*(numproc-1));
  }
#endif

  progname = basename(*argv);
  if(argc == 1){
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
  verbose = 0; wavelet = 0; pad = 0; fixed = 0;
  pairid = 0; n2id = 0;
  parinfo.ntrig = 331;
  parinfo.picture_threshold = DEFAULT_PICTHRESH; 
  parinfo.boundary_threshold = DEFAULT_BNDTHRESH;
  parinfo.efactor = DEFAULT_EFACTOR;
  strcpy(path_to_db,DEFAULT_PATH);
  strcpy(parinfo.source_name,"UNKNOWN");
  parinfo.source_ra = 0.0; parinfo.source_dec = 0.0;
  parinfo.xo = 0.0; parinfo.yo = 0.0;
/*************************************************************************/
/* get command line options                                              */
/*************************************************************************/
  getoptions(argc, argv);
  filename = *(argv+optind);

  if(myid == master){ /* master does IO */
/*************************************************************************/
/* open input file                                                       */
/*************************************************************************/
    infile_id = Hopen(filename,DFACC_RDONLY,0);
    if(infile_id == -1){
      fprintf(stderr,"%s: can\'t open %s\n",progname,filename);
      ierr = 1; /* set error flag */
    }else{
      Vstart(infile_id);

/*************************************************************************/
/* set runinfo fields to read                                            */
/*************************************************************************/
      vruninfo_ref = VSfind(infile_id,"10M Run Information");
      vruninfo_id = VSattach(infile_id,vruninfo_ref,"r");
      VSsetfields(vruninfo_id,
		  "VERSION,NADC,UTC_START,SOURCE_NAME,SOURCE_RA,SOURCE_DEC");
      runinfobuf = (uchar8 *)malloc(2*sizeof(int32)+3*sizeof(float64)+
				    80*sizeof(char8));

/*************************************************************************/
/* read run information                                                  */
/*************************************************************************/
      printf("%s: reading %s ...\n",progname,filename);
      VSread(vruninfo_id,runinfobuf,1,FULL_INTERLACE);
      memcpy(&(parinfo.version),runinfobuf,4);
      memcpy(&(parinfo.nadc),runinfobuf+4,4);
      memcpy(&(parinfo.utc_start),runinfobuf+8,8);
      if(!strcmp(parinfo.source_name,"UNKNOWN")){
	memcpy(parinfo.source_name,runinfobuf+16,80);
	parinfo.source_name[24] = '\0'; 
      }
      if(parinfo.source_ra == 0.0)
	memcpy(&(parinfo.source_ra),runinfobuf+96,8);
      if(parinfo.source_dec == 0.0)
	memcpy(&(parinfo.source_dec),runinfobuf+104,8);
      parinfo.npmt = get_npmt(parinfo.nadc);

/*************************************************************************/
/* set event fields to read                                              */
/*************************************************************************/
      vevent_ref = VSfind(infile_id,"10M Event");
      vevent_id = VSattach(infile_id,vevent_ref,"r");
      VSsetfields(vevent_id,
		  "CODE,GPSTIME,OSCTIME,LIVETIME,ADC,ELEVATION,AZIMUTH");

/*************************************************************************/
/* set param info and output to screen                                   */
/*************************************************************************/
      strncpy(fileid,basename(filename),8); fileid[8] = '\0';
/*************************************************************************/
/* idate yymmdd                                                          */
/*************************************************************************/
      sla_djcl_(&(parinfo.utc_start),&iyy,&imm,&idd,&dfr,&j);
      sla_calyd_(&iyy,&imm,&idd,&nyy,&ndd,&j);
      if(iyy < 2000)
	parinfo.idate = (iyy-1900)*10000+imm*100+idd;
      else
	parinfo.idate = (iyy-2000)*10000+imm*100+idd;
      printf("   run date: %06d\n",parinfo.idate);
/*************************************************************************/
/* itime hhmmss (utc)                                                    */
/*************************************************************************/
      j = 1; sla_dd2tf_(&j,&dfr,&cr,itime,1);
      parinfo.itime = itime[0]*10000 + itime[1]*100 + itime[2];
      printf("   run time: %2dh %2dm %2ds\n",itime[0],itime[1],itime[2]);
/*************************************************************************/
/* sidereal time                                                         */
/*************************************************************************/
      siderealtime = sla_gmst_(&(parinfo.utc_start)) - WHIPLONG;
      siderealtime = sla_dranrm_(&siderealtime);
      j = 1; sla_dr2tf_(&j,&siderealtime,&cr,ra,1);
      printf("   sidereal time: %c%dh %dm %d.%ds\n",cr,
	     ra[0],ra[1],ra[2],ra[3]);
/*************************************************************************/
/* source name                                                           */
/*************************************************************************/
      printf("   source: %s\n",parinfo.source_name);
/*************************************************************************/
/* source coordinates current epoch                                      */
/*************************************************************************/
      j = 1; sla_dr2tf_(&j,&(parinfo.source_ra),&cr,ra,1);
      sla_dr2af_(&j,&(parinfo.source_dec),&cd,dec,1);
      epoch = (double)iyy + (double)(ndd/365.0);
      printf("   ra,dec: %c%dh %dm %d.%ds,%c%dd %dm %d.%ds (J%.1lf)\n",
	     cr,ra[0],ra[1],ra[2],ra[3],cd,dec[0],dec[1],dec[2],dec[3],epoch);
/*************************************************************************/
/* source coordinates J2000                                              */
/*************************************************************************/
      ra2000 = parinfo.source_ra; dec2000 = parinfo.source_dec; j2000 = 2000.0;
      sla_preces_("FK4",&epoch,&j2000,&(ra2000),&(dec2000),3);
      j = 1; sla_dr2tf_(&j,&ra2000,&cr,ra,1);
      sla_dr2af_(&j,&dec2000,&cd,dec,1);
      printf("   ra,dec: %c%dh %dm %d.%ds,%c%dd %dm %d.%ds (J2000.0)\n",
	     cr,ra[0],ra[1],ra[2],ra[3],cd,dec[0],dec[1],dec[2],dec[3]);
/*************************************************************************/
/* camera info                                                           */
/*************************************************************************/
      printf("   no. channels: %d\n",parinfo.nadc);
      printf("   no. pmts:     %d\n",parinfo.npmt);
      printf("   no. trigger:  %d\n",parinfo.ntrig);
    }
  }

#ifdef MPI
  MPI_Bcast(&parinfo,1,parinfo_datatype,master,MPI_COMM_WORLD);
  MPI_Bcast(&ierr,1,MPI_INT,master,MPI_COMM_WORLD);
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
  adc = (short *)malloc(parinfo.nadc*sizeof(short));
  peds = (float *)malloc(sizeof(float)*parinfo.nadc);
  pedvars = (float *)malloc(sizeof(float)*parinfo.nadc);
  pairpeds = (float *)malloc(sizeof(float)*parinfo.nadc);
  pairpedvars = (float *)malloc(sizeof(float)*parinfo.nadc);
  maxpedvars = (float *)malloc(sizeof(float)*parinfo.nadc);
  gains = (float *)malloc(sizeof(float)*parinfo.nadc);
  picture = (int *)malloc(sizeof(int)*parinfo.npmt);
  boundary = (int *)malloc(sizeof(int)*parinfo.npmt);
  parinfo.tubesoff = (int *)malloc(sizeof(int)*parinfo.npmt);
  pairtubesoff = (int *)malloc(sizeof(int)*parinfo.npmt);
  fadc = (float *)malloc(sizeof(float)*parinfo.npmt);
  for(j=0;j < parinfo.npmt;j++){
    parinfo.tubesoff[j] = 0;
    pairtubesoff[j] = 0;
  }

/*************************************************************************/
/* pmt coordinates and neighbours                                        */
/*************************************************************************/
  pmt = cwhipplecoords(parinfo.npmt);
  pmt_neighbours = cwhippleneighbours(parinfo.npmt);
  pmt_spacing = 7.0;/*(float)sqrt((double)(pmt.wxdeg(1)*pmt.wxdeg(1) + 
		      pmt.wydeg(1)*pmt.wydeg(1)))*60.0;*/

  if(myid == master){
/*************************************************************************/
/* retrieve database info                                                */
/*************************************************************************/
    printf("   reading database information from %s\n",path_to_db);
/*************************************************************************/
/* pedestals                                                             */
/*************************************************************************/
    if(!getpeds(fileid,parinfo.idate,peds,pedvars,path_to_db)){
      printf("%s: can\'t find pedestals for %s %06d\n",progname,fileid,
	     parinfo.idate);
      ierr = 1;
    }
    if(!pairid && pad){
      if(!getpair(fileid,parinfo.idate,parinfo.pairid,path_to_db)){
	printf("%s: can\'t find padding pair id for %s\n",progname,fileid);
	ierr = 1;
      }
    }
    if(pad){
      if(!getpeds(parinfo.pairid,parinfo.idate,pairpeds,pairpedvars,
		  path_to_db)){
	printf("%s: can\'t find pedestals for %s %06d\n",progname,
	       parinfo.pairid,parinfo.idate);
	ierr = 1;
      }
    }
/*************************************************************************/
/* nitrogen gains                                                        */
/*************************************************************************/
    if(!n2id){
      if(!getnlist(fileid,parinfo.idate,parinfo.n2id,path_to_db)){
	printf("%s: can\'t find nitrogen id for %s\n",progname,fileid);
	ierr = 1;
      }
    }
    if(!getn2gains(parinfo.n2id,parinfo.idate,gains,path_to_db)){
      printf("%s: can\'t find nitrogen gain factors for %s %06d\n",
	     progname,parinfo.n2id,parinfo.idate);
      ierr = 1;
    }
/*************************************************************************/
/* tubes off                                                             */
/*************************************************************************/
    getoff(fileid,parinfo.idate,parinfo.tubesoff,path_to_db);
    for(j = 0; (parinfo.tubesoff[j] != 0) && (j < parinfo.npmt);j++);
    parinfo.ntubesoff = j;
    if(pad){
      getoff(parinfo.pairid,parinfo.idate,pairtubesoff,path_to_db);
      for(i = 0; (pairtubesoff[i] != 0) && (i < parinfo.npmt);i++){
	for(k = 0;k < parinfo.ntubesoff;k++){
	  if(parinfo.tubesoff[k] == pairtubesoff[i])
	    break;
	  if(k == (parinfo.ntubesoff-1))
	    parinfo.tubesoff[j++] = pairtubesoff[i];
	}
      }
      parinfo.ntubesoff = j;
    }
  }

#ifdef MPI
  MPI_Bcast(&ierr,1,MPI_INT,master,MPI_COMM_WORLD);
  MPI_Bcast(peds,parinfo.nadc,MPI_FLOAT,master,MPI_COMM_WORLD);
  MPI_Bcast(pedvars,parinfo.nadc,MPI_FLOAT,master,MPI_COMM_WORLD);
  if(pad){
    MPI_Bcast(pairpeds,parinfo.nadc,MPI_FLOAT,master,MPI_COMM_WORLD);
    MPI_Bcast(pairpedvars,parinfo.nadc,MPI_FLOAT,master,MPI_COMM_WORLD);
  }
  MPI_Bcast(gains,parinfo.nadc,MPI_FLOAT,master,MPI_COMM_WORLD);
  MPI_Bcast(parinfo.tubesoff,parinfo.npmt,MPI_INT,master,MPI_COMM_WORLD);
  MPI_Bcast(&parinfo.ntubesoff,1,MPI_INT,master,MPI_COMM_WORLD);
#endif
  if(ierr){ /* quit on error */
#ifdef MPI
    MPI_Finalize();
#endif
    return 1;
  }

  if(myid == master){
/*************************************************************************/
/* parameterization info                                                 */
/*************************************************************************/
    if(pad)
      strcpy(parinfo.padlabel,"ON");
    else
      strcpy(parinfo.padlabel,"OFF");
    printf("   padding is %3s\n",parinfo.padlabel);
    if(pad)
      printf("   padding pair is %8s\n",parinfo.pairid);
    strcpy(parinfo.gainlabel,"N2");
    printf("   nitrogen gain id is %8s\n",parinfo.n2id);
    printf("   picture/boundary thresholds = %5.2f/%5.2f\n",
	   parinfo.picture_threshold,parinfo.boundary_threshold);
    printf("   elongation factor = %5.2f\n",parinfo.efactor);

/*************************************************************************/
/* open outfile                                                          */
/*************************************************************************/
    strcpy(HDF_outname,fileid); HDF_outname[8]='\0';
    strncat(HDF_outname,"p.hdf",5);
    printf("   output file is %s\n",HDF_outname);
    outfile_id = Hopen(HDF_outname,DFACC_CREATE,0);
    if(outfile_id == -1){
      fprintf(stderr,"%s: can\'t open %s\n",progname,HDF_outname);
      ierr = 1;
    }else{
      Vstart(outfile_id);

      vparinfo_id = VSattach(outfile_id, -1, "w");
      VSfdefine(vparinfo_id,vparinfo_name[0],DFNT_INT32,   1); /*SW VERSION */
      VSfdefine(vparinfo_id,vparinfo_name[1],DFNT_INT32,   1); /*NADC       */
      VSfdefine(vparinfo_id,vparinfo_name[2],DFNT_FLOAT64, 1); /*UTC_START  */
      VSfdefine(vparinfo_id,vparinfo_name[3],DFNT_INT32,   1); /*IDATE      */
      VSfdefine(vparinfo_id,vparinfo_name[4],DFNT_INT32,   1); /*ITIME      */
      VSfdefine(vparinfo_id,vparinfo_name[5],DFNT_CHAR8,  80); /*SOURCENAME */
      VSfdefine(vparinfo_id,vparinfo_name[6],DFNT_FLOAT64, 1); /*SOURCE RA  */
      VSfdefine(vparinfo_id,vparinfo_name[7],DFNT_FLOAT64, 1); /*SOURCE DEC */
      VSfdefine(vparinfo_id,vparinfo_name[8],DFNT_CHAR8,   3); /*GAIN LABEL */
      VSfdefine(vparinfo_id,vparinfo_name[9],DFNT_CHAR8,   9); /*N2ID       */
      VSfdefine(vparinfo_id,vparinfo_name[10],DFNT_CHAR8,  4); /*PAD LABEL  */
      VSfdefine(vparinfo_id,vparinfo_name[11],DFNT_INT32,  1); /*NTUBESOFF  */
      ntubesoff = parinfo.ntubesoff > 0 ? parinfo.ntubesoff : 1;
      VSfdefine(vparinfo_id,vparinfo_name[12],DFNT_INT32,ntubesoff);
      VSfdefine(vparinfo_id,vparinfo_name[13],DFNT_FLOAT32,1); /*PICTURE    */
      VSfdefine(vparinfo_id,vparinfo_name[14],DFNT_FLOAT32,1); /*BOUNDARY   */
      VSfdefine(vparinfo_id,vparinfo_name[15],DFNT_FLOAT32,1); /*EFACTOR    */
      VSsetname(vparinfo_id,"10M Parameter Info");
      VSsetclass(vparinfo_id,"Parameter Info");
      VSsetfields(vparinfo_id,vparinfo_name_str);
      parinfobuf = (uint8 *)malloc((5+ntubesoff)*sizeof(int32) +
				   3*sizeof(float64) + 3*sizeof(float32) +
				   96*sizeof(char8));

      vpar_id = VSattach(outfile_id, -1, "w");
      VSfdefine(vpar_id,vpar_name[0],DFNT_FLOAT64, 1); /*GPSTIME */
      VSfdefine(vpar_id,vpar_name[1],DFNT_FLOAT64, 1); /*OSCTIME */
      VSfdefine(vpar_id,vpar_name[2],DFNT_FLOAT64, 1); /*LIVETIME*/
      VSfdefine(vpar_id,vpar_name[3],DFNT_FLOAT32, 1); /*LENGTH  */
      VSfdefine(vpar_id,vpar_name[4],DFNT_FLOAT32, 1); /*WIDTH   */
      VSfdefine(vpar_id,vpar_name[5],DFNT_FLOAT32, 1); /*MISS    */
      VSfdefine(vpar_id,vpar_name[6],DFNT_FLOAT32, 1); /*DISTANCE*/
      VSfdefine(vpar_id,vpar_name[7],DFNT_FLOAT32, 1); /*AZWIDTH */
      VSfdefine(vpar_id,vpar_name[8],DFNT_FLOAT32, 3); /*FRAC    */
      VSfdefine(vpar_id,vpar_name[9],DFNT_INT32,   1); /*SIZE    */
      VSfdefine(vpar_id,vpar_name[10],DFNT_INT32,  3); /*LOC     */
      VSfdefine(vpar_id,vpar_name[11],DFNT_FLOAT32,3); /*MAX     */
      VSfdefine(vpar_id,vpar_name[12],DFNT_FLOAT32,1); /*ASYMM   */
      VSfdefine(vpar_id,vpar_name[13],DFNT_FLOAT32,1); /*XC      */
      VSfdefine(vpar_id,vpar_name[14],DFNT_FLOAT32,1); /*YC      */
      VSfdefine(vpar_id,vpar_name[15],DFNT_FLOAT32,1); /*XO      */
      VSfdefine(vpar_id,vpar_name[16],DFNT_FLOAT32,1); /*YO      */
      VSsetname(vpar_id,"10M Parameters");
      VSsetclass(vpar_id,"Parameters");
      VSsetfields(vpar_id,vpar_name_str);

/*************************************************************************/
/* write param info                                                      */
/*************************************************************************/
      memcpy(parinfobuf,&(parinfo.version),4);
      memcpy(parinfobuf+4,&(parinfo.nadc),4);
      memcpy(parinfobuf+8,&(parinfo.utc_start),8);
      memcpy(parinfobuf+16,&(parinfo.idate),4);
      memcpy(parinfobuf+20,&(parinfo.itime),4);
      memcpy(parinfobuf+24,parinfo.source_name,80);
      memcpy(parinfobuf+104,&(parinfo.source_ra),8);
      memcpy(parinfobuf+112,&(parinfo.source_dec),8);
      memcpy(parinfobuf+120,parinfo.gainlabel,3);
      memcpy(parinfobuf+123,parinfo.n2id,9);
      memcpy(parinfobuf+132,parinfo.padlabel,4);
      memcpy(parinfobuf+136,&(parinfo.ntubesoff),4);
      memcpy(parinfobuf+140,parinfo.tubesoff,ntubesoff*4);
      memcpy(parinfobuf+140+ntubesoff*4,&(parinfo.picture_threshold),4);
      memcpy(parinfobuf+140+ntubesoff*4+4,&(parinfo.boundary_threshold),4);
      memcpy(parinfobuf+140+ntubesoff*4+8,&(parinfo.efactor),4);
      VSwrite(vparinfo_id, (uchar8 *)parinfobuf, 1, FULL_INTERLACE);

      printf("   ok ... processing   ");
      fflush(stdout);
    }
  }

#ifdef MPI
  MPI_Bcast(&ierr,1,MPI_INT,master,MPI_COMM_WORLD);
#endif
  if(ierr){ /* quit on error */
#ifdef MPI
    MPI_Finalize();
#endif
    return 1;
  }
  parbuf = (uint8 **)malloc((numproc-1)*sizeof(uint8 *));
  for(i=0;i < numproc;i++){
    parbuf[i] = (uint8 *)malloc(3*sizeof(float64) + 16*sizeof(float32) +
				4*sizeof(int32));
  }
  eventbuf = (uchar8 *)malloc((parinfo.nadc+3)*sizeof(int16) + 
			      3*sizeof(float64));

/*************************************************************************/
/* MAIN LOOP                                                             */
/*************************************************************************/
  i = 1, k = 0; l = 1;
  while(1){
/*************************************************************************/
/* master does IO                                                        */
/*************************************************************************/
    if(myid == master){
      if(VSread(vevent_id,eventbuf,1,FULL_INTERLACE)==1){ /* read event */
	memcpy(&code,eventbuf,2);
	if(code & 4){
/*************************************************************************/
/* if no MPI or only 1 processor no need to send data as master is slave */
/*************************************************************************/
#ifdef MPI	  
	  if(numproc > 1){
/*************************************************************************/
/* initially give each processor an event to work on                     */
/*************************************************************************/
	    if(i < numproc){
	      MPI_Send(eventbuf,(parinfo.nadc+3)*2+24,MPI_BYTE,i,1,
		       MPI_COMM_WORLD);
	      MPI_Irecv(parbuf[i-1],104,MPI_BYTE,i,MPI_ANY_TAG,
			MPI_COMM_WORLD,&requests[i-1]);
	      i++;
	    }else{
/*************************************************************************/
/* give receive event to whoever is finished                             */
/*************************************************************************/
	      MPI_Waitany(numproc-1,requests,&proc,&status[0]);
	      proc++;
	      if(status[0].MPI_TAG){
		k++;
		VSwrite(vpar_id, (uchar8 *)parbuf[proc-1],1, FULL_INTERLACE);
		if(wavelet){
		  switch(k%8){
		  case 0: case 4:
		    printf("\b|");
		    break;
		  case 1: case 5:
		    printf("\b/");
		    break;
		  case 2: case 6:
		    printf("\b-");
		    break;
		  case 3: case 7:
		    printf("\b\\");
		  }
		  fflush(stdout);
		}
	      }
	      MPI_Send(eventbuf,(parinfo.nadc+3)*2+24,MPI_BYTE,proc,1,
		       MPI_COMM_WORLD);
	      MPI_Irecv(parbuf[proc-1],104,MPI_BYTE,proc,MPI_ANY_TAG,
			MPI_COMM_WORLD,&requests[proc-1]);
	    }
	  }
	  //JPF There was a misplaced bracket here... it needed to be moved
	  //outside the #endif below
#endif
      }}else{
/*************************************************************************/
/* EOF - wait for all processors to return their results then send stop  */
/*************************************************************************/
#ifdef MPI
	if(numproc > 1)
	  MPI_Waitall(numproc-1,requests,status);
	for(i=1;i < numproc;i++){
	  if(status[i-1].MPI_TAG)
	    VSwrite(vpar_id, (uchar8 *)parbuf[i-1], 1, FULL_INTERLACE);
	  MPI_Send(eventbuf,0,MPI_BYTE,i,0,MPI_COMM_WORLD); /* stop slaves */
	}
	MPI_Finalize();
#endif
	printf("\n");
/*************************************************************************/
/* finished parameterizing - close files                                 */
/*************************************************************************/
	VSdetach(vruninfo_id);
	VSdetach(vevent_id);
	Vend(infile_id);
	Hclose(infile_id);

	VSdetach(vparinfo_id);
	VSdetach(vpar_id);
	Vend(outfile_id);
	Hclose(outfile_id);

	printf("done ...\n");

	return 0;
      }
    }
/*************************************************************************/
/* slave                                                                 */
/*************************************************************************/
    if(myid != 0 || (numproc == 1 && code & 4)){
#ifdef MPI
/*************************************************************************/
/* wait for an event to work on                                          */
/*************************************************************************/
      if(numproc > 1){
	MPI_Recv(eventbuf,(parinfo.nadc+3)*2+24,MPI_BYTE,master,
		 MPI_ANY_TAG,MPI_COMM_WORLD,&status[0]);
	if(!status[0].MPI_TAG){ /* received stop */
	  MPI_Finalize();
	  return 0;
	}
      }
#endif
/*************************************************************************/
/* prepare data for parameterizing                                       */
/*************************************************************************/
      memcpy(&(par.gpstime),eventbuf+2,8);
      memcpy(&(par.osctime),eventbuf+10,8);
      memcpy(&(par.livetime),eventbuf+18,8);
      memcpy(adc,eventbuf+26,2*parinfo.nadc);
      memcpy(&elevation,eventbuf+26+2*parinfo.nadc,2);
      memcpy(&azimuth,eventbuf+26+2*parinfo.nadc+2,2);

      for(j=0;j < parinfo.npmt;j++)
	fadc[j] = (float)(abs(adc[j]));

      cimpad(fadc,pedvars,pairpedvars,maxpedvars);
      if(wavelet){
	sub_macro_wv_analysis__(&(parinfo.npmt),fadc,
			      &pmt_spacing,
			      peds,maxpedvars,
			      pmt.wxdeg,pmt.wydeg,
			      &j,&l);
	l = 0;
      }else{
	cimclean(fadc,peds,maxpedvars,picture,boundary);
      }
      cimgain(fadc,gains);

      if(verbose){
	for(j=0;j < parinfo.npmt;j++){
	  printf("%4d %4d %6.2f %6.2f %6.2f %6.2f\n",
		 j+1,adc[j],peds[j],pedvars[j],
		 gains[j],fadc[j]);
	}
	/* printf("Time of event: %lf\n",par.gpstime); */
	/* getchar(); */
      }

/*************************************************************************/
/* parameterize                                                          */
/*************************************************************************/
      if(get_hillas_parameters(fadc,&par,parinfo.xo,parinfo.yo,0)){
	memcpy(parbuf[0],&(par.gpstime),8);
	memcpy(parbuf[0]+8,&(par.osctime),8);
	memcpy(parbuf[0]+16,&(par.livetime),8);
	memcpy(parbuf[0]+24,&(par.length),4);
	memcpy(parbuf[0]+28,&(par.width),4);
	memcpy(parbuf[0]+32,&(par.miss),4);
	memcpy(parbuf[0]+36,&(par.distance),4);
	memcpy(parbuf[0]+40,&(par.azwidth),4);
	memcpy(parbuf[0]+44,par.frac,12);
	memcpy(parbuf[0]+56,&(par.size),4);
	memcpy(parbuf[0]+60,par.location,12);
	memcpy(parbuf[0]+72,par.maximum,12);
	memcpy(parbuf[0]+84,&(par.asymm),4);
	memcpy(parbuf[0]+88,&(par.xc),4);
	memcpy(parbuf[0]+92,&(par.yc),4);
	memcpy(parbuf[0]+96,&(par.xo),4);
	memcpy(parbuf[0]+100,&(par.yo),4);
/*************************************************************************/
/* MPI - send back results or write if only one processor                */
/*************************************************************************/
#ifdef MPI
	if(numproc > 1){
	  MPI_Send(parbuf[0],104,MPI_BYTE,master,1,MPI_COMM_WORLD);
	}else{
	  VSwrite(vpar_id, (uchar8 *)parbuf[0], 1, FULL_INTERLACE);
	}
#else
/*************************************************************************/
/* no MPI do write results                                               */
/*************************************************************************/
	VSwrite(vpar_id, (uchar8 *)parbuf[0], 1, FULL_INTERLACE);
#endif
      }else{
/*************************************************************************/
/* parameterize failed send failed TAG                                   */
/*************************************************************************/
#ifdef MPI
	if(numproc > 1)
	  MPI_Send(parbuf[0],0,MPI_BYTE,master,0,MPI_COMM_WORLD);
#else
	;
#endif
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
  int c, j;
  int rahours, ramins, decdegs, decmins;
  double rasecs, decsecs;

  while(1){
/*************************************************/
/* Note: getopt is POSIX.1 compliant thus should */
/*       available on all newer workstations.    */
/*       getopt_long would be better as it sorts */
/*       argv into options -> non-options, hence */
/*       the options could appear anywhere on    */
/*       the command line                        */
/*************************************************/  
    c = getopt(argc, argv, "b:c:d:e:fg:hi:n:pr:s:t:vw");
    switch(c){
    case 'b':
      if(sscanf(optarg,"%f,%f",&(parinfo.picture_threshold),
		&(parinfo.boundary_threshold)) != 2)
	printf("%s: invalid option\n",progname);
      break;
    case 'c':
      if(sscanf(optarg,"%f,%f",&(parinfo.xo),&(parinfo.yo)) != 2)
	printf("%s: invalid option\n",progname);
      break;
    case 'd':
      if(sscanf(optarg,"%lf",&(parinfo.source_dec)) != 1)
	printf("%s: invalid option\n",progname);
      decdegs = (int)(parinfo.source_dec/10000.0);
      decmins = (int)((parinfo.source_dec - (double)decdegs*10000.0)/100.0);
      decsecs = parinfo.source_dec - (double)decdegs*10000.0 - 
	(double)decmins*100.0;
      sla_daf2r_(&decdegs,&decmins,&decsecs,&(parinfo.source_dec),&j);
      break;
    case 'e':
      if(sscanf(optarg,"%f",&(parinfo.efactor)) != 1)
	printf("%s: invalid option\n",progname);
      break;
    case 'f':
      fixed = 1;
      break;
    case 'g':
      strcpy(path_to_db,optarg);
      break;
    case 'h':
      usage();
#ifdef MPI
      MPI_Finalize();
#endif
      exit(EXIT_SUCCESS);
    case 'i':
      if(sscanf(optarg,"%s",parinfo.pairid) != 1)
	printf("%s: invalid option\n",progname);
      pairid = 1;
      break;
    case 'n':
      if(sscanf(optarg,"%s",parinfo.n2id) != 1)
	printf("%s: invalid option\n",progname);
      n2id = 1;
      break;
    case 'p':
      pad = 1;
      break;
    case 'r':
      if(sscanf(optarg,"%lf",&(parinfo.source_ra)) != 1)
	printf("%s: invalid option\n",progname);
      rahours = (int)(parinfo.source_ra/10000.0);
      ramins = (int)((parinfo.source_ra - (double)rahours*10000.0)/100.0);
      rasecs = parinfo.source_ra - (double)rahours*10000.0 -
	(double)ramins*100.0;
      sla_dtf2r_(&rahours,&ramins,&rasecs,&(parinfo.source_ra),&j);
      break;
    case 's':
      if(sscanf(optarg,"%s",parinfo.source_name) != 1)
	printf("%s: invalid option\n",progname);
      break;
    case 't':
      if(sscanf(optarg,"%d",&(parinfo.ntrig)) != 1)
	printf("%s: invalid option\n",progname);
      break;
    case 'v':
      verbose = 1;
      break;
    case 'w':
      wavelet = 1;
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
  printf("   -b %%f,%%f  picture,Boundary threshold [%.2f,%.2f]\n",
	 DEFAULT_PICTHRESH,DEFAULT_BNDTHRESH);
  printf("   -d %%f  Declination of source [ddmmss.s]\n");
  printf("   -e %%f  Elongation factor [%.2f]\n",
	 DEFAULT_EFACTOR);
  printf("   -f     Fixed two dimensional analysis, no derotation\n");
  printf("   -g %%s  path to database files [%s]\n",DEFAULT_PATH);
  printf("   -h     Help\n");
  printf("   -i %%s  padding pair Id\n");
  printf("   -n %%s  Nitrogen id\n");
  printf("   -p     software Padding\n");
  printf("   -r %%f  Right ascension of source [hhmmss.s]\n");
  printf("   -s %%s  Source name\n");
  printf("   -t %%d  number of Trigger tubes [%d]\n",parinfo.ntrig);
  printf("   -v     Verbose text output of all events\n");
  printf("   -w     clean images using Wavelets\n");
  return;
}

/*************************************************************************/
/* function: impad                                                       */
/* passed:   see below                                                   */
/* returns:  0 - fail, 1 - ok                                            */
/* purpose:  pads images                                                 */
/*************************************************************************/
int cimpad(float *event, float *pedvars, float *pairpedvars, float *maxpedvars)
{
  static int idum = -13;
  int k;

/*************************************************************************/
/* software padding                                                      */
/*************************************************************************/
  if(pad){
    for(k=0;k < parinfo.npmt;k++){
      maxpedvars[k] = (pairpedvars[k] > pedvars[k]) ? 
	pairpedvars[k] : pedvars[k]; 
      if(pairpedvars[k] > pedvars[k]) /* add some noise */
	event[k] += (float)(sqrt((double)(pairpedvars[k]*pairpedvars[k] - 
				  pedvars[k]*pedvars[k])))*gasdev(&idum);
    }
  }else{
    for(k=0;k < parinfo.npmt;k++)
      maxpedvars[k] = pedvars[k];
  }

  return 1;
}

/*************************************************************************/
/* function: imclean                                                     */
/* passed:   see below                                                   */
/* returns:  0 - fail, 1 - ok                                            */
/* purpose:  cleans images                                               */
/*************************************************************************/
int cimclean(float *event, float *peds, float *maxpedvars, 
	     int *picture, int *boundary)
{
  int k, i;

  for(k=0;k < parinfo.npmt;k++) /* Subtract peds */
    event[k] = event[k] - peds[k];
  k = 0;
  while(parinfo.tubesoff[k] != 0) /* turn off tubes */
    event[parinfo.tubesoff[k++]-1] = 0.0;
  for(k=0;k < parinfo.npmt;k++) /* picture tubes */
    picture[k] = (event[k] > parinfo.picture_threshold*maxpedvars[k]) ? 1 : 0;
  for(k=0;k < parinfo.npmt;k++){ /* boundary tubes */
    boundary[k] = 0;
    if(event[k] > parinfo.boundary_threshold*maxpedvars[k]){
      for(i=0;i < NNEIGHBOURS;i++)
	  if((pmt_neighbours[k][i]>0)&&
	     picture[pmt_neighbours[k][i]-1]) boundary[k] = 1;
    }
  }
  for(k=0;k < parinfo.npmt;k++){
    if(!picture[k]&&!boundary[k]) event[k] = 0.0; /*zero tubes not in picture*/
  }
  return 1;
}

/*************************************************************************/
/* function: cimgain                                                     */
/* passed:   adc values and gain factors                                 */
/* returns:  0 - fail, 1 - ok                                            */
/* purpose:  applys gain normalizations                                  */
/*************************************************************************/
int cimgain(float *event, float *gains)
{
  int k;

  for(k=0;k < parinfo.npmt;k++){
    event[k] *= gains[k]; /* apply gains */
  }
  return 1;
}

/*************************************************************************/
/* function: get_hillas_parameters                                       */
/* passed:   see below                                                   */
/* returns:  0 - fail, 1 - ok                                            */
/* purpose:  parameterizes images                                        */
/*************************************************************************/
int get_hillas_parameters(float *event, struct par_s *par, float xo, float yo, 
			  int offsetflag)
{
  float sumsig, sumxsig, sumx2sig, sumysig, sumy2sig, sumxysig;
  float sumx3sig, sumx2ysig, sumxy2sig, sumy3sig;
  float xmean, x2mean, ymean, y2mean, xymean, xmean2, ymean2, meanxy;
  float x3mean, x2ymean, xy2mean, y3mean;
  float sdevx2, sdevy2, sdevxy, wxbyev, wybyev;
  float sdevx3, sdevx2y, sdevxy2, sdevy3;
  float d, a, delta, z, u, v, psi, cpsi, spsi;
  float slope, yint;
  double qa, qb, qc;
  float temp, xrel, yrel, xt, yt;
  int i;
  struct par_s temppar;

  sumsig = 0.0; sumxsig = 0.0; sumx2sig = 0.0; 
  sumysig = 0.0; sumy2sig = 0.0; sumxysig = 0.0;
  sumx3sig = 0.0; sumx2ysig = 0.0;
  sumxy2sig = 0.0; sumy3sig = 0.0;

  for(i=0;i < parinfo.npmt;i++){
    if(event[i] > 0.0){
      /* Allows parameters to be calculated relative to coordinates xo,yo */
      xrel = pmt.wxdeg[i] - xo;
      yrel = pmt.wydeg[i] - yo;
      wxbyev = xrel*event[i];
      wybyev = yrel*event[i];
      sumsig += event[i];
      sumxsig += wxbyev;
      sumx2sig += xrel*wxbyev;
      sumysig += wybyev;
      sumy2sig += yrel*wybyev;
      sumxysig += xrel*wybyev;
      sumx3sig += xrel*xrel*wxbyev;
      sumx2ysig += xrel*xrel*wybyev;
      sumxy2sig += xrel*yrel*wybyev;
      sumy3sig += yrel*yrel*wybyev;
    }
  }
  if(sumsig <= 0.0) return 0;

  xmean = sumxsig/sumsig;
  x2mean = sumx2sig/sumsig;
  ymean = sumysig/sumsig;
  y2mean = sumy2sig/sumsig;
  xymean = sumxysig/sumsig;
  x3mean = sumx3sig/sumsig;
  x2ymean = sumx2ysig/sumsig;
  xy2mean = sumxy2sig/sumsig;
  y3mean = sumy3sig/sumsig;
  xmean2 = xmean*xmean;
  ymean2 = ymean*ymean;
  meanxy = xmean*ymean;

  sdevx2 = x2mean - xmean2;
  sdevy2 = y2mean - ymean2;
  sdevxy = xymean - meanxy;
  sdevx3 = x3mean - 3.0*xmean*x2mean + 2.0*xmean2*xmean;
  sdevy3 = y3mean - 3.0*ymean*y2mean + 2.0*ymean2*ymean;
  sdevx2y = x2ymean - x2mean*ymean - 2.0*xymean*xmean + 2.0*xmean2*ymean;
  sdevxy2 = xy2mean - y2mean*xmean - 2.0*xymean*ymean + 2.0*xmean*ymean2;

  d = sdevy2 - sdevx2;
  temp = d*d + 4.0*sdevxy*sdevxy;
  z = (temp <= 0.0) ? 0.0 : sqrt(temp);

/******************************* Distance ***********************************/
  temp = xmean2 + ymean2;
  par->distance = (temp <= 0.0) ? 0.0 : sqrt(temp);
/****************************************************************************/
/******************************* Length *************************************/
  temp = (sdevx2+sdevy2 + z)/2.0;
  par->length = (temp <= 0.0) ? 0.0 : sqrt(temp);
/****************************************************************************/
/******************************* Width **************************************/
  temp = (sdevy2+sdevx2 - z)/2.0;
  par->width = (temp <= 0.0) ? 0.0 : sqrt(temp);
/****************************************************************************/
/******************************* Miss ***************************************/
  if(z == 0.0){
    par->miss = par->distance;
  }else{
    u = 1.0 + d/z; v = 2.0 - u;
    temp = (u*xmean2 + v*ymean2)/2.0 - meanxy*(2.0*sdevxy/z);
    par->miss = (temp <= 0.0) ? 0.0 : sqrt(temp);
  }
/****************************************************************************/
/******************** equation of straight line y = mx + b ******************/
  if(!offsetflag){
    temp = d*d+4.0*sdevxy*sdevxy;
    if(sdevxy == 0.0)
	slope = 999.;   /* vertical line essentially */
    else
      slope = (temp <= 0.0) ? 0.0 : (d+sqrt(temp))/2.0/sdevxy;
    yint = ymean - slope*xmean;
  }
/****************************************************************************/
/******************************* Asymmetry **********************************/
  if(par->length == 0.0){
    par->asymm = 0.0;
  }else{
    psi = atan2((d+z)*ymean + 2.0*sdevxy*xmean,2.0*sdevxy*ymean - (d-z)*xmean);
    cpsi=cos(psi);
    spsi=sin(psi);
    par->asymm=(sdevx3*cpsi*cpsi*cpsi+3.0*sdevx2y*spsi*cpsi*cpsi+
	   3.0*sdevxy2*cpsi*spsi*spsi+sdevy3*spsi*spsi*spsi);
    par->asymm=(par->asymm < 0.0) ? 
      -1.0*exp(log(fabs(par->asymm))/3.0)/par->length : 
      exp(log(fabs(par->asymm))/3.0)/par->length;
  }
/****************************************************************************/
/******************************* Centroid ***********************************/
  par->xc = xmean; 
  par->yc = ymean;
/****************************************************************************/
  d = y2mean - x2mean;
  temp = d*d + 4.0*xymean*xymean;
  z = (temp <= 0.0) ? 0.0 : sqrt(temp);
/******************************* Azwidth ************************************/
  temp = (x2mean + y2mean - z)/2.0;
  par->azwidth = (temp  <= 0.0) ? 0.0 : sqrt(temp);
/****************************************************************************/
  getmax(event,parinfo.npmt,par->maximum,par->location);
  getfrac(par->maximum,par->frac,sumsig);
  getmax(event,parinfo.ntrig,par->maximum,par->location);
  par->size = (int)sumsig;
  if(!offsetflag){
/*****************************************************************************
     Point of origin section - for 2 d analysis
     distance along line from centroid to point of origin 
*****************************************************************************/
    if(par->length == 0.00)
     a = 999.;
    else
      a = parinfo.efactor - parinfo.efactor*par->width/par->length;
    qa = (double)slope*(double)slope + 1.0;
    qb = 2.0*(double)slope*(double)yint - 
      2.0*(double)ymean*(double)slope - 2.0*(double)xmean;
    qc = (double)xmean*(double)xmean + (double)ymean*(double)ymean - 
      2.0*(double)ymean*(double)yint + 
      (double)yint*(double)yint - (double)a*(double)a;
    xt = (float)((-qb-sqrt(qb*qb - 4.0*qa*qc))/2.0/qa);
    yt = slope*xt + yint;
    /* Use asymmetry to decide which side of the ellipse to use for 
       point of origin */
    get_hillas_parameters(event,&temppar,xt,yt,1);
    if(temppar.asymm > 0.0){
      par->xo = xt;
      par->yo = yt;
    }else{
      par->xo= (-qb+sqrt(qb*qb - 4.0*qa*qc))/2.0/qa;
      par->yo = slope*par->xo + yint;
    }
/*****************************************************************************
  The above positions are in camera coordinates. In order to combine
  data from different runs (i.e. different siderial times) they must
  be de-rotated to the same hour angle.
*****************************************************************************/
    if(!fixed){
      derot(&(par->xo),&(par->yo),par->gpstime);
    }
  }

  return 1;
}

/*************************************************************************/
/* function: derot_angle                                                 */
/* passed:   nothing                                                     */
/* returns:  angle to derot coordinates (uses actual telescope el and az)*/
/* purpose:  coordinate derotation                                       */
/*************************************************************************/
double derot_angle(double gpstime)
{
  double altitude, azimuth;
  double hourangle, derotangle;
  static double latitude = WHIPLAT;

  hourangle = sla_gmst_(&gpstime) - WHIPLONG;
  hourangle = sla_dranrm_(&hourangle) - parinfo.source_ra;
  
  sla_de2h_(&hourangle,&(parinfo.source_dec),&latitude,&azimuth,&altitude);
  derotangle = atan2(-1.0*COSLATITUDE*sin(azimuth),
		     (cos(altitude)*SINLATITUDE - sin(altitude)*cos(azimuth)));

  return derotangle;
}

/*************************************************************************/
/* function: derot                                                       */
/* passed:   x, y coords to derotate                                     */
/* returns:  nothing                                                     */
/* purpose:  does actual derotation                                      */
/*************************************************************************/
void derot(float *x, float *y, double gpstime)
{
  double tmp_x, tmp_y, theta;

  theta = derot_angle(gpstime);
  tmp_x = (double)(*x)*cos(theta)-(double)(*y)*sin(theta);
  tmp_y = (double)(*x)*sin(theta)+(double)(*y)*cos(theta);

  *x = (float)tmp_x;
  *y = (float)tmp_y;

  return;
}

/*************************************************************************/
/* function: getfrac                                                     */
/* passed:   adcs and size                                               */
/* returns:  frac1-3                                                     */
/* purpose:  finds above                                                 */
/*************************************************************************/
void getfrac(float *x, float *frac, float total)
{
  int i;

  frac[0] = x[0]/total;
  for(i=1;i < 3;i++)
    frac[i] = frac[i-1] + x[i]/total;
    
}
    
/*************************************************************************/
/* function: getmax                                                      */
/* passed:   adc values and number of adcs                               */
/* returns:  max1-3 loc1-3                                               */
/* purpose:  finds above                                                 */
/*************************************************************************/
void getmax(float *x, int N, float max[], int loc[])
{
  int i;

  max[0] = max[1] = max[2] = 0.0;
  loc[0] = loc[1] = loc[2] = 0;

  for(i = 0;i < N;i++){
    if(x[i] > max[2]){
      if(x[i] > max[1]){
	if(x[i] > max[0]){
	  max[2] = max[1];
	  max[1] = max[0];
	  max[0] = x[i];
	  loc[2] = loc[1];
	  loc[1] = loc[0];
	  loc[0] = i;
	}else{
	  max[2] = max[1];
	  max[1] = x[i];
	  loc[2] = loc[1];
	  loc[1] = i;
	}
      }else{
	max[2] = x[i];
	loc[2] = i;
      }
    }
  }

  return;
}
/* Thats all folks */
