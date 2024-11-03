/****************************************************************************/
/* mcut: performs cuts on parameterised files                               */
/* Written by R.W. Lessard 960606                                           */
/* Usage: mcut [options] files                                              */
/*        where files are parameterised                                     */
/*        files. See usage for option list.                                 */
/* Patchlist:                                                               */
/*            Y2K compliant 12/27/1999                                      */
/*            user definable output 01/28/2000                              */
/*            MPI 04/11/00                                                  */
/****************************************************************************/
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
#include "mbook.h"
#include "hdf.h"
#ifdef MPI
#include "mpi.h"
#endif

/*************************************************************************/
/* defines                                                               */
/*************************************************************************/
#define MAX_ALPHA               90.0
#define MAX_THETA               2.0
#define MAX_FOV                 1.95
#define RAW                     0
#define TRIG                    1
#define SHAPE                   2
#define ALPHA                   3
#define APERTURE                4

/* Default cut values */
#define DEFAULT_UPPER_WIDTH     0.12
#define DEFAULT_LOWER_WIDTH     0.05
#define DEFAULT_UPPER_LENGTH    0.25
#define DEFAULT_LOWER_LENGTH    0.12
#define DEFAULT_UPPER_DISTANCE  1.00
#define DEFAULT_LOWER_DISTANCE  0.40
#define DEFAULT_UPPER_ALPHA     15.0
#define DEFAULT_LOWER_ASYMM     -999.9
#define DEFAULT_UPPER_THETA     0.22
#define DEFAULT_FRAC            0.975
#define DEFAULT_TRIGGER_0       30
#define DEFAULT_TRIGGER_1       30
#define DEFAULT_TRIGGER_2       0
#define DEFAULT_LOWER_SIZE      0
#define DEFAULT_UPPER_SIZE      999999
#define DEFAULT_LENGTH_SIZE     0.00040
#define DEFAULT_DURATION        999999.0
#define DEFAULT_ALPHA_BIN       5.0
#define DEFAULT_NPHASE_BINS     50
#define DEFAULT_THETA_BIN       0.02
#define DEFAULT_TWOD_BIN        0.1
#define DEFAULT_IMAGE_XC        0.0
#define DEFAULT_IMAGE_YC        0.0

/*************************************************/
/* storage for the parameters                    */
/*************************************************/
struct par_s {
  double gpstime;
  double solar_bary;
  double orbit_bary;
  double phase;
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
} par;

/*************************************************/
/* storage for the photon data                   */
/*************************************************/
struct photon_s {
  double time[100000];
  double solar_bary[100000];
  double orbit_bary[100000];
  double phase[100000];
  double energy[100000];
};

/*************************************************/
/* storage for the results                       */
/*************************************************/
struct output_results_s {
  long passed[5];
  MHist alpha;
  MHist theta;
  MHist raw_image;
  MHist smoothed_image;
  struct photon_s photon;
  MHist phase;
} output_results;

/*************************************************/
/* storage for super cuts                        */
/*************************************************/
struct cut_s {
  int lower_size;
  int upper_size;
  int trigger[3];
  float upper_width;
  float lower_width;
  float upper_length;
  float lower_length;
  float upper_distance;
  float lower_distance;
  float upper_alpha;
  float lower_asymm;
  float upper_theta;
  float frac;
  float length_size;
  double duration;
  int nphase_bins;
  float alpha_bin;
  float theta_bin;
  float twod_bin;
  float image_xc;
  float image_yc;
} cut = {DEFAULT_LOWER_SIZE,DEFAULT_UPPER_SIZE,
	 DEFAULT_TRIGGER_0,DEFAULT_TRIGGER_1,DEFAULT_TRIGGER_2,
	 DEFAULT_UPPER_WIDTH,DEFAULT_LOWER_WIDTH,
	 DEFAULT_UPPER_LENGTH,DEFAULT_LOWER_LENGTH,
	 DEFAULT_UPPER_DISTANCE,DEFAULT_LOWER_DISTANCE,
	 DEFAULT_UPPER_ALPHA,DEFAULT_LOWER_ASYMM,
	 DEFAULT_UPPER_THETA,DEFAULT_FRAC,
	 DEFAULT_LENGTH_SIZE,DEFAULT_DURATION,DEFAULT_NPHASE_BINS,
	 DEFAULT_ALPHA_BIN,DEFAULT_THETA_BIN,DEFAULT_TWOD_BIN,
	 DEFAULT_IMAGE_XC,DEFAULT_IMAGE_YC};

/*************************************************/
/* storage for periodic analysis                   */
/*************************************************/
struct pulsar_parameters_s{
  double freq;
  double freqdot;
  double epoch;
} pulsar;

struct binary_parameters_s{
  double Po;
  double Po_dot;
  double ecc;
  double asini;
  double Omega;
  double Omega_dot;
  double gamma;
  double to;
} binary;

/*************************************************/
/* storage for the header                        */
/*************************************************/
struct parinfo_s{
  char source_name[80];
  double utc_start;
  double source_ra;
  double source_dec;
  double duration;
  int idate;
  int itime;
} parinfo;

char *progname;
char *filename;
int myid = 0, master = 0;

/*************************************************/
/* Options                                       */
/*************************************************/
extern int optind;
int extended_cuts;
int small_cuts;
int gaussian;
int phase;
int clock;
int orbit;
int selector;

/*************************************************/
/* Function declarations                         */
/*************************************************/
void cut_parameters();
void cut_parameters_s();
void cut_parameters_e();
void getoptions();
void usage();
void get_phase(double, MHist *);
double get_orbit(double);
double bary(double, int, int, double, double);
double mfun2g(double, double);
double mfun2b(double, double);

int main(int argc, char *argv[])
{
  int i, j, k, l;
  
  int32 infile_id, vparinfo_ref, vparinfo_id, vpar_ref, vpar_id;
  uchar8 *parinfobuf, *parbuf;

  struct output_results_s output_tmp;

  int numproc = 1, ierr = 0, proc = 1;

#ifdef MPI
  /* parameter structure */
  int par_blockcounts[3] = {5,16,4};
  MPI_Datatype par_types[3] = {MPI_DOUBLE,MPI_FLOAT,MPI_INT};
  MPI_Aint par_displs[3];
  MPI_Datatype par_datatype;

  /* info structure */
  int parinfo_blockcounts[3] = {80,4,2};
  MPI_Datatype parinfo_types[3] = {MPI_CHAR,MPI_DOUBLE,MPI_INT};
  MPI_Aint parinfo_displs[3];
  MPI_Datatype parinfo_datatype;

  MPI_Status stat;

/*************************************************************************/
/* initialize MPI                                                        */
/*************************************************************************/
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numproc);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);

/*************************************************************************/
/* define some MPI data types                                            */
/*************************************************************************/
  MPI_Address(&par.gpstime,&par_displs[0]);
  MPI_Address(&par.length,&par_displs[1]);
  MPI_Address(&par.size,&par_displs[2]);
  for(i=2;i >= 0; i--)
    par_displs[i] -= par_displs[0];
  MPI_Type_struct(3,par_blockcounts,par_displs,par_types,&par_datatype);
  MPI_Type_commit(&par_datatype);

  MPI_Address(&parinfo.source_name,&parinfo_displs[0]);
  MPI_Address(&parinfo.utc_start,&parinfo_displs[1]);
  MPI_Address(&parinfo.idate,&parinfo_displs[2]);
  for(i=2;i >= 0; i--)
    parinfo_displs[i] -= parinfo_displs[0];
  MPI_Type_struct(3,parinfo_blockcounts,parinfo_displs,
		  parinfo_types,&parinfo_datatype);
  MPI_Type_commit(&parinfo_datatype);
#endif

  progname = (char *)basename(*argv);
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
  gaussian = 0; phase = 0; clock = 1; orbit = 0; 
  extended_cuts = 0; small_cuts = 0; selector = ALPHA;


/*************************************************************************/
/* get command line options                                              */
/*************************************************************************/
  getoptions(argc,argv);

  filename = *(argv+optind);
  mbook1(&(output_results.phase),cut.nphase_bins,0.0,2.0);
  mbook1(&(output_results.alpha),
	 (int)(MAX_ALPHA/cut.alpha_bin+0.1),0.0,MAX_ALPHA);
  mbook1(&(output_results.theta),
	 (int)(MAX_THETA/cut.theta_bin+0.1),0.0,MAX_THETA);
  mbook2(&(output_results.raw_image),
	 (int)(2*MAX_FOV/cut.twod_bin+0.1),-MAX_FOV,MAX_FOV,
	 (int)(2*MAX_FOV/cut.twod_bin+0.1),-MAX_FOV,MAX_FOV);
  mbook2(&(output_results.smoothed_image),
	 (int)(2*MAX_FOV/cut.twod_bin+0.1),-MAX_FOV,MAX_FOV,
	 (int)(2*MAX_FOV/cut.twod_bin+0.1),-MAX_FOV,MAX_FOV);

  if((myid == master) && (numproc > 1)){
    mbook1(&(output_tmp.phase),cut.nphase_bins,0.0,2.0);
    mbook1(&(output_tmp.alpha),
	   (int)(MAX_ALPHA/cut.alpha_bin+0.1),0.0,MAX_ALPHA);
    mbook1(&(output_tmp.theta),
	   (int)(MAX_THETA/cut.theta_bin+0.1),0.0,MAX_THETA);
    mbook2(&(output_tmp.raw_image),
	   (int)(2*MAX_FOV/cut.twod_bin+0.1),-MAX_FOV,MAX_FOV,
	   (int)(2*MAX_FOV/cut.twod_bin+0.1),-MAX_FOV,MAX_FOV);
    mbook2(&(output_tmp.smoothed_image),
	   (int)(2*MAX_FOV/cut.twod_bin+0.1),-MAX_FOV,MAX_FOV,
	   (int)(2*MAX_FOV/cut.twod_bin+0.1),-MAX_FOV,MAX_FOV);
  }
  output_results.passed[0] = 0; output_results.passed[1] = 0;
  output_results.passed[2] = 0; output_results.passed[3] = 0;
  output_results.passed[4] = 0;

  if(myid == master){ /* master does IO */
/*************************************************************************/
/* open HDF file                                                         */
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
      vparinfo_ref = VSfind(infile_id,"10M Parameter Info");
      vparinfo_id = VSattach(infile_id,vparinfo_ref,"r");
      VSsetfields(vparinfo_id,"UTC_START,IDATE,ITIME,SOURCE_NAME,SOURCE_RA,SOURCE_DEC");
      parinfobuf = (uint8 *)malloc(2*sizeof(int32) + 3*sizeof(float64) +
				   80*sizeof(char8));

/*************************************************************************/
/* read param info                                                       */
/*************************************************************************/
      VSread(vparinfo_id,parinfobuf,1,FULL_INTERLACE);
      memcpy(&(parinfo.utc_start),parinfobuf,8);
      memcpy(&(parinfo.idate),parinfobuf+8,4);
      memcpy(&(parinfo.itime),parinfobuf+12,4);
      memcpy(parinfo.source_name,parinfobuf+16,80);
      if(phase == 0){
	memcpy(&(parinfo.source_ra),parinfobuf+96,8);
	memcpy(&(parinfo.source_dec),parinfobuf+104,8);
      }
      parinfo.duration = 0.0;

/*************************************************************************/
/* set parameter fields to read                                          */
/*************************************************************************/
      vpar_ref = VSfind(infile_id,"10M Parameters");
      vpar_id = VSattach(infile_id,vpar_ref,"r");
      VSsetfields(vpar_id,"GPSTIME,OSCTIME,LIVETIME,LENGTH,WIDTH,MISS,DISTANCE,AZWIDTH,FRAC,SIZE,LOC,MAX,ASYMM,XC,YC,XO,YO");
      parbuf = (uint8 *)malloc(3*sizeof(float64) + 16*sizeof(float32) +
			       4*sizeof(int32));
    }
  }

#ifdef MPI
  MPI_Bcast(&parinfo,1,parinfo_datatype,master,MPI_COMM_WORLD);
  MPI_Bcast(&ierr,1,MPI_INT,master,MPI_COMM_WORLD);    /* send ierr to all */
#endif
  if(ierr){ /* quit on error */
#ifdef MPI
    MPI_Finalize();
#endif
    return 1;
  }

  if(myid == 0)
    printf("%s: processing %s ...\n",progname,filename);

/*************************************************************************/
/* MAIN LOOP                                                             */
/*************************************************************************/
  while(1){
/*************************************************************************/
/* master                                                                */
/*************************************************************************/
    if(myid == master){ /* master reads file and sends data to slaves */
      if((VSread(vpar_id,parbuf,1,FULL_INTERLACE) == 1)
	 && (parinfo.duration < cut.duration)){
	memcpy(&(par.gpstime),parbuf,8);
	memcpy(&(par.osctime),parbuf+8,8);
	memcpy(&(par.livetime),parbuf+16,8);
	memcpy(&(par.length),parbuf+24,4);
	memcpy(&(par.width),parbuf+28,4);
	memcpy(&(par.miss),parbuf+32,4);
	memcpy(&(par.distance),parbuf+36,4);
	memcpy(&(par.azwidth),parbuf+40,4);
	memcpy(par.frac,parbuf+44,12);
	memcpy(&(par.size),parbuf+56,4);
	memcpy(par.location,parbuf+60,12);
	memcpy(par.maximum,parbuf+72,12);
	memcpy(&(par.asymm),parbuf+84,4);
	memcpy(&(par.xc),parbuf+88,4);
	memcpy(&(par.yc),parbuf+92,4);
	memcpy(&(par.xo),parbuf+96,4);
	memcpy(&(par.yo),parbuf+100,4);

	if(parinfo.duration < par.livetime/60.0)
	  parinfo.duration = par.livetime/60.0;
#ifdef MPI	  
	if(numproc > 1){
	  MPI_Send(&par,1,par_datatype,proc,1,MPI_COMM_WORLD);
	  if(proc < numproc-1) proc++; else proc = 1;
	}
#endif
      }else{
#ifdef MPI	
/*************************************************************************/
/* collect data from slaves                                              */
/*************************************************************************/
	k = 0, l = 0;
	for(j=1;j < numproc;j++){
	  MPI_Send(&par,0,par_datatype,j,0,MPI_COMM_WORLD);
	  /* events passed cuts */
	  MPI_Recv(output_tmp.passed,5,MPI_LONG,j,1,MPI_COMM_WORLD,&stat);
	  /* photon data */
	  MPI_Recv(output_tmp.photon.time,output_tmp.passed[selector],
		   MPI_DOUBLE,j,1,MPI_COMM_WORLD,&stat);
	  MPI_Recv(output_tmp.photon.energy,output_tmp.passed[selector],
		   MPI_DOUBLE,j,1,MPI_COMM_WORLD,&stat);

	  /* distributions */
	  MPI_Recv(output_tmp.alpha.matrixptr,output_results.alpha.nibins,
		   MPI_DOUBLE,j,1,MPI_COMM_WORLD,&stat);
	  MPI_Recv(output_tmp.theta.matrixptr,output_results.theta.nibins,
		   MPI_DOUBLE,j,1,MPI_COMM_WORLD,&stat);
	  MPI_Recv(output_tmp.raw_image.matrixptr,
		   output_results.raw_image.nibins*
		   output_results.raw_image.njbins,
		   MPI_DOUBLE,j,1,MPI_COMM_WORLD,&stat);
	  MPI_Recv(output_tmp.smoothed_image.matrixptr,
		   output_results.smoothed_image.nibins*
		   output_results.smoothed_image.njbins,
		   MPI_DOUBLE,j,1,MPI_COMM_WORLD,&stat);

	  /* periodic analysis results */
	  if(phase){
	    MPI_Recv(output_tmp.photon.solar_bary,output_tmp.passed[selector],
		     MPI_DOUBLE,j,1,MPI_COMM_WORLD,&stat);
	    MPI_Recv(output_tmp.photon.orbit_bary,output_tmp.passed[selector],
		     MPI_DOUBLE,j,1,MPI_COMM_WORLD,&stat);
	    MPI_Recv(output_tmp.photon.phase,output_tmp.passed[selector],
		     MPI_DOUBLE,j,1,MPI_COMM_WORLD,&stat);

	    MPI_Recv(output_tmp.phase.matrixptr,
		     output_results.phase.nibins,
		     MPI_DOUBLE,j,1,MPI_COMM_WORLD,&stat);
	  }
/*************************************************************************/
/* accumulate results                                                    */
/*************************************************************************/
	  /* events passed cuts */
	  for(i=0;i < 5;i++)
	    output_results.passed[i] += output_tmp.passed[i];

	  /* photon data */
	  for(i=0;i < output_tmp.passed[selector];i++){
	    output_results.photon.time[k] = 
	      output_tmp.photon.time[i];
	    output_results.photon.energy[k++] = 
	      output_tmp.photon.energy[i];
	    if(phase){
	      output_results.photon.solar_bary[k-1] = 
		output_tmp.photon.solar_bary[i];
	      output_results.photon.orbit_bary[k-1] = 
		output_tmp.photon.orbit_bary[i];
	      output_results.photon.phase[k-1] = 
		output_tmp.photon.phase[i];
	    }
	  }

	  /* distributions */
	  for(i=0;i < output_results.alpha.nibins;i++)
	    *(output_results.alpha.matrixptr+i) += 
	      *(output_tmp.alpha.matrixptr+i);
	  for(i=0;i < output_results.theta.nibins;i++)
	    *(output_results.theta.matrixptr+i) += 
	      *(output_tmp.theta.matrixptr+i);
	  for(i=0;i < output_results.raw_image.nibins*
		output_results.raw_image.njbins;i++)
	    *(output_results.raw_image.matrixptr+i) += 
	      *(output_tmp.raw_image.matrixptr+i);
	  for(i=0;i < output_results.smoothed_image.nibins*
		output_results.smoothed_image.njbins;i++)
	    *(output_results.smoothed_image.matrixptr+i) += 
	      *(output_tmp.smoothed_image.matrixptr+i);
	  if(phase){
	    for(i=0;i < output_results.phase.nibins;i++)
	      *(output_results.phase.matrixptr+i) += 
		*(output_tmp.phase.matrixptr+i);
	  }
	}
#endif
	VSdetach(vparinfo_id);
	VSdetach(vpar_id);
	Vend(infile_id);
	Hclose(infile_id);

#ifdef HAVE_MATLAB
	MATLAB_output();
#endif
#ifdef HAVE_USER_MCUT_OUTPUT
	USER_MCUT_output();
#endif
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
    if(myid != 0 || numproc == 1){
#ifdef MPI
      if(numproc > 1){
	MPI_Recv(&par,1,par_datatype,master,
		 MPI_ANY_TAG,MPI_COMM_WORLD,&stat);
	if(!stat.MPI_TAG){ /* received stop, send results */
/*************************************************************************/
/* send results to master                                                */
/*************************************************************************/
	  /* events passed cuts */
	  MPI_Send(output_results.passed,5,MPI_LONG,master,1,MPI_COMM_WORLD);

	  /* photon data */
	  MPI_Send(output_results.photon.time,output_results.passed[selector],
		   MPI_DOUBLE,master,1,MPI_COMM_WORLD);
	  MPI_Send(output_results.photon.energy,output_results.passed[selector],
		   MPI_DOUBLE,master,1,MPI_COMM_WORLD);

	  /* distributions */
	  MPI_Send(output_results.alpha.matrixptr,
		   output_results.alpha.nibins,
		   MPI_DOUBLE,master,1,MPI_COMM_WORLD);
	  MPI_Send(output_results.theta.matrixptr,
		   output_results.theta.nibins,
		   MPI_DOUBLE,master,1,MPI_COMM_WORLD);
	  MPI_Send(output_results.raw_image.matrixptr,
		   output_results.raw_image.nibins*
		   output_results.raw_image.njbins,
		   MPI_DOUBLE,master,1,MPI_COMM_WORLD);
	  MPI_Send(output_results.smoothed_image.matrixptr,
		   output_results.smoothed_image.nibins*
		   output_results.smoothed_image.njbins,
		   MPI_DOUBLE,master,1,MPI_COMM_WORLD);

          /* periodic analysis results */
	  if(phase){
	    MPI_Send(output_results.photon.solar_bary,
		     output_results.passed[selector],
		     MPI_DOUBLE,master,1,MPI_COMM_WORLD);
	    MPI_Send(output_results.photon.orbit_bary,
		     output_results.passed[selector],
		     MPI_DOUBLE,master,1,MPI_COMM_WORLD);
	    MPI_Send(output_results.photon.phase,
		     output_results.passed[selector],
		     MPI_DOUBLE,master,1,MPI_COMM_WORLD);

	    MPI_Send(output_results.phase.matrixptr,
		     output_results.phase.nibins,
		     MPI_DOUBLE,master,1,MPI_COMM_WORLD);
	  }
	  MPI_Finalize();
	  return 0;
	}
      }
#endif
      if(extended_cuts)
	cut_parameters_e();
      else
	cut_parameters();
    }
  }
}    

/*************************************************/
/* supercuts                                     */
/*************************************************/  
void cut_parameters()
{
  double alpha, theta, time;
  
/*************************************************/
/* increment raw event counter                   */
/*************************************************/  
  output_results.passed[RAW]++;
/*************************************************/
/* size, trigger, length/size and cosmic ray cuts*/
/*************************************************/
  if((par.size >= cut.lower_size) && (par.size <= cut.upper_size) &&
     (par.maximum[0] >= cut.trigger[0]) &&
     (par.maximum[1] >= cut.trigger[1]) &&
     (par.maximum[2] >= cut.trigger[2]) &&
     ((par.length/(float)par.size) <= cut.length_size) &&
     (par.frac[2] < cut.frac)){
    output_results.passed[TRIG]++;
    alpha = asin(par.miss/par.distance)*180.0/M_PI;
    theta = sqrt((par.xo-cut.image_xc)*(par.xo-cut.image_xc) + 
		(par.yo-cut.image_yc)*(par.yo-cut.image_yc));
    if((par.length < cut.upper_length) && 
       (par.length > cut.lower_length) &&
       (par.width < cut.upper_width) && 
       (par.width > cut.lower_width)){
/*************************************************/
/* passed shape cut                              */
/*************************************************/  
      output_results.passed[SHAPE]++;
/*************************************************/
/* alpha                                         */
/*************************************************/  
      if((par.distance < cut.upper_distance) && 
	 (par.distance > cut.lower_distance)){
/*************************************************/
/* passed distance cut                           */
/*************************************************/  
	if(par.asymm > cut.lower_asymm){
	  mfill(&(output_results.alpha),alpha,alpha);
	  if(alpha < cut.upper_alpha){
/*************************************************/
/* passed shape and alpha cut                    */
/*************************************************/  
	    output_results.passed[ALPHA]++;
	    if(selector == ALPHA){ /* photon data selected by ALPHA cut */
	      if(clock == 1)
		time = par.gpstime;
	      else
		time = (double)((int)(par.gpstime)) + par.osctime/86400.0L;
	      output_results.photon.time[output_results.passed[ALPHA]-1] 
		= time;
	      output_results.photon.energy[output_results.passed[ALPHA]-1] 
		= (double)par.size;
	      if(phase){
		get_phase(time,&output_results.phase);
		output_results.photon.solar_bary[output_results.passed[ALPHA]-1] = par.solar_bary;
		output_results.photon.orbit_bary[output_results.passed[ALPHA]-1] =	par.orbit_bary;
		output_results.photon.phase[output_results.passed[ALPHA]-1] = par.phase;
	      }
	    }
	  }
	}
      }
/*************************************************/
/* aperture                                      */
/*************************************************/  
      mfill(&(output_results.theta),theta,theta);
      mfill(&(output_results.raw_image),(double)(par.xo),(double)(par.yo));
      if(gaussian)
	mfun(&(output_results.smoothed_image),mfun2g);
      else
	mfun(&(output_results.smoothed_image),mfun2b);
      if(theta < cut.upper_theta){
/*************************************************/
/* passed shape and theta cut                   */
/*************************************************/  
	output_results.passed[APERTURE]++;
	if(selector == APERTURE){
	  if(clock == 1)
	    time = par.gpstime;
	  else
	    time = (double)((int)(par.gpstime)) + par.osctime/86400.0L;
	  output_results.photon.time[output_results.passed[APERTURE]-1] = time;
	  output_results.photon.energy[output_results.passed[APERTURE]-1] = 
	    (double)par.size;
	  if(phase){
	    get_phase(time,&output_results.phase);
	    output_results.photon.solar_bary[output_results.passed[APERTURE]-1] = par.solar_bary;
	    output_results.photon.orbit_bary[output_results.passed[APERTURE]-1] = par.orbit_bary;
	    output_results.photon.phase[output_results.passed[APERTURE]-1] = par.phase;
	  }
	}
      }
    }
  }else{
    if(small_cuts)
      cut_parameters_s();
  }

  return;
}

/*************************************************/
/* small cuts - hardwired                        */
/*************************************************/
void cut_parameters_s()
{
  double alpha, theta, time;

/*************************************************/
/* size, trigger, length/size and cosmic ray cuts*/
/*************************************************/
  if((par.size >= 0) &&
     (par.maximum[0] >= 30) &&
     (par.maximum[1] >= 30) &&
     (par.maximum[2] >= 0) &&
     ((par.length/(float)par.size) <= 0.00042) &&
     (par.frac[2] < cut.frac)){
    output_results.passed[TRIG]++;
    alpha = asin(par.miss/par.distance)*180.0/M_PI;
    theta = sqrt((par.xo-cut.image_xc)*(par.xo-cut.image_xc) + 
		(par.yo-cut.image_yc)*(par.yo-cut.image_yc));
    if((par.length < 0.21) && 
       (par.length > 0.12) &&
       (par.width < 0.12) && 
       (par.width > 0.03)){
/*************************************************/
/* passed shape cut                              */
/*************************************************/  
      output_results.passed[SHAPE]++;
/*************************************************/
/* alpha                                         */
/*************************************************/  
      if((par.distance < 1.00) && 
	 (par.distance > 0.36)){
/*************************************************/
/* passed distance cut                           */
/*************************************************/  
	if(par.asymm > cut.lower_asymm){
	  mfill(&(output_results.alpha),alpha,alpha);
	  if(alpha < 16.0){
/*************************************************/
/* passed shape and alpha cut                    */
/*************************************************/  
	    output_results.passed[ALPHA]++;
	    if(selector == ALPHA){
	      if(clock == 1)
		time = par.gpstime;
	      else
		time = (double)((int)(par.gpstime)) + par.osctime/86400.0L;
	      output_results.photon.time[output_results.passed[ALPHA]-1] = 
		time;
	      output_results.photon.energy[output_results.passed[ALPHA]-1] = (double)par.size;
	      if(phase){
		get_phase(time,&output_results.phase);
		output_results.photon.solar_bary[output_results.passed[ALPHA]-1] =	par.solar_bary;
		output_results.photon.orbit_bary[output_results.passed[ALPHA]-1] =	par.orbit_bary;
		output_results.photon.phase[output_results.passed[ALPHA]-1] = par.phase;
	      }
	    }
	  }
	}
      }

/*************************************************/
/* aperture                                      */
/*************************************************/  
      mfill(&(output_results.theta),theta,theta);
      mfill(&(output_results.raw_image),(double)(par.xo),(double)(par.yo));
      if(gaussian)
	mfun(&(output_results.smoothed_image),mfun2g);
      else
	mfun(&(output_results.smoothed_image),mfun2b);
      if(theta < cut.upper_theta){
/*************************************************/
/* passed shape and theta cut                   */
/*************************************************/  
	output_results.passed[APERTURE]++;
	if(selector == APERTURE){
	  if(clock == 1)
	    time = par.gpstime;
	  else
	    time = (double)((int)(par.gpstime)) + par.osctime/86400.0L;
	  output_results.photon.time[output_results.passed[APERTURE]-1] = time;
	  output_results.photon.energy[output_results.passed[APERTURE]-1] = (double)par.size;
	  if(phase){
	    get_phase(time,&output_results.phase);
	    output_results.photon.solar_bary[output_results.passed[APERTURE]-1] = par.solar_bary;
	    output_results.photon.orbit_bary[output_results.passed[APERTURE]-1] = par.orbit_bary;
	    output_results.photon.phase[output_results.passed[APERTURE]-1] = par.phase;
	  }
	}
      }
    }
  }

  return;
}

/*************************************************/
/* extended-supercuts                            */
/*************************************************/  
void cut_parameters_e()
{
  double LC, WC;
  double alpha, theta, time;

/*************************************************/
/* increment raw event counter                   */
/*************************************************/  
  output_results.passed[RAW]++;
/*************************************************/
/* size, trigger and cosmic ray cuts             */
/*************************************************/
  if((par.size >= cut.lower_size) && (par.size <= cut.upper_size) &&
     (par.maximum[0] >= cut.trigger[0]) &&
     (par.maximum[1] >= cut.trigger[1]) &&
     (par.maximum[2] >= cut.trigger[2]) &&
     ((par.length/(float)par.size) <= cut.length_size) &&
     (par.frac[2] < cut.frac)){
    output_results.passed[TRIG]++;
    alpha = asin(par.miss/par.distance)*180.0/M_PI;
    theta = sqrt((par.xo-cut.image_xc)*(par.xo-cut.image_xc) + 
		(par.yo-cut.image_yc)*(par.yo-cut.image_yc));
    LC = fabs((double)(par.length) - 0.144 - 0.020*log((double)par.size));
    WC = fabs((double)(par.width) + 0.022 - 0.023*log((double)par.size));
    if((LC < 0.068) && (WC < 0.048)){
/*************************************************/
/* passed shape cut                              */
/*************************************************/  
      output_results.passed[SHAPE]++;
/*************************************************/
/* alpha                                         */
/*************************************************/  
      if((par.distance < 1.0) && 
	 (par.distance > 0.6)){
/*************************************************/
/* passed distance cut                           */
/*************************************************/  
	mfill(&(output_results.alpha),alpha,alpha);
	if((alpha < (13.5 - 0.558*log((double)par.size) + 9.6))){
/*************************************************/
/* passed shape and alpha cut                    */
/*************************************************/  
	  output_results.passed[ALPHA]++;
	  if(selector == ALPHA){
	    if(clock == 1)
	      time = par.gpstime;
	    else
	      time = (double)((int)(par.gpstime)) + par.osctime/86400.0L;
	    output_results.photon.time[output_results.passed[ALPHA]-1] = time;
	    output_results.photon.energy[output_results.passed[ALPHA]-1] = (double)par.size;
	    if(phase){
	      get_phase(time,&output_results.phase);
	      output_results.photon.solar_bary[output_results.passed[ALPHA]-1] = par.solar_bary;
	      output_results.photon.orbit_bary[output_results.passed[ALPHA]-1] = par.orbit_bary;
	      output_results.photon.phase[output_results.passed[ALPHA]-1] = par.phase;
	    }
	  }
	}
      }
/*************************************************/
/* aperture                                      */
/*************************************************/  
      mfill(&(output_results.theta),theta,theta);
      mfill(&(output_results.raw_image),(double)(par.xo),(double)(par.yo));
      if(gaussian)
	mfun(&(output_results.smoothed_image),mfun2g);
      else
	mfun(&(output_results.smoothed_image),mfun2b);
      if(theta < cut.upper_theta){
/*************************************************/
/* passed shape and theta cut                   */
/*************************************************/  
	output_results.passed[APERTURE]++;
	if(selector == APERTURE){
	  if(clock == 1)
	    time = par.gpstime;
	  else
	    time = (double)((int)(par.gpstime)) + par.osctime/86400.0L;
	  output_results.photon.time[output_results.passed[APERTURE]-1] = time;
	  output_results.photon.energy[output_results.passed[APERTURE]-1] = (double)par.size;
	  if(phase){
	    get_phase(time,&output_results.phase);
	    output_results.photon.solar_bary[output_results.passed[APERTURE]-1] = par.solar_bary;
	    output_results.photon.orbit_bary[output_results.passed[APERTURE]-1] = par.orbit_bary;
	    output_results.photon.phase[output_results.passed[APERTURE]-1] = par.phase;
	  }
	}
      }
    }
  }

  return;
}

/*************************************************/
/* command line options                          */
/*************************************************/  
void getoptions(int argc, char **argv)
{
  extern char *optarg;
  int c;
  float tmp;

  while(1){
/*************************************************/
/* Note: getopt is POSIX.1 compliant thus should */
/*       available on all newer workstations.    */
/*       getopt_long would be better as it sorts */
/*       argv into options -> non-options, hence */
/*       the options could appear anywhere on    */
/*       the command line                        */
/*************************************************/  
    c = getopt(argc, argv, "a:b:c:d:ef:ghk:l:mn:o:r:p:s:t:u:y:w:z:");
    switch(c){
    case 'a':
      if(sscanf(optarg,"%f",&cut.upper_alpha) != 1)
	printf("%s: invalid alpha option\n",progname);
      break;
    case 'b':
      if(sscanf(optarg,"%f,%f,%f",&cut.alpha_bin,&cut.theta_bin,
		&cut.twod_bin) != 3)
	printf("%s: invalid bin option\n",progname);
      break;
    case 'c':
      if(sscanf(optarg,"%f,%f",&cut.image_xc,&cut.image_yc) != 2)
	printf("%s: invalid aperture center option\n",progname);
      break;
    case 'd':
      if(sscanf(optarg,"%f,%f",&cut.lower_distance,&cut.upper_distance) == 2){
	if(cut.upper_distance < cut.lower_distance){
	  tmp = cut.upper_distance;
	  cut.upper_distance = cut.lower_distance;
	  cut.lower_distance = tmp;
	}
      } else printf("%s: invalid distance range option\n",progname);
      break;
    case 'e':
      extended_cuts = 1;
      break;
    case 'f':
      if(sscanf(optarg,"%f",&cut.frac) != 1)
	printf("%s: invalid duration option\n",progname);
      break;
    case 'g':
      gaussian = 1;
      break;
    case 'h':
      usage();
#ifdef MPI
      MPI_Finalize();
#endif
      exit(EXIT_SUCCESS);
    case 'k':
      if(sscanf(optarg,"%lf",&cut.duration) != 1)
	printf("%s: invalid duration option\n",progname);
      break;
    case 'l':
      if(sscanf(optarg,"%f,%f",&cut.lower_length,&cut.upper_length) == 2){
	if(cut.upper_length < cut.lower_length){
	  tmp = cut.upper_length;
	  cut.upper_length = cut.lower_length;
	  cut.lower_length = tmp;
	}
      } else printf("%s: invalid length range option\n",progname);
      break;
    case 'm':
      small_cuts = 1;
      break;
    case 'n':
      if(sscanf(optarg,"%d",&cut.nphase_bins) != 1)
	printf("%s: invalid phase bin option\n",progname);
      break;
    case 'o':
      if(sscanf(optarg,"%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf",
		&binary.Po,&binary.Po_dot,&binary.ecc,&binary.asini,
		&binary.Omega,&binary.Omega_dot,&binary.gamma,&binary.to) != 8)
	printf("%s: invalid binary option\n",progname);
      orbit = 1;
      break;
    case 'p':
      if(sscanf(optarg,"%lf,%lf,%lf,%lf,%lf,%d",
		&parinfo.source_ra,&parinfo.source_dec,
		&pulsar.freq,&pulsar.freqdot,&pulsar.epoch,
		&clock) != 6)
	printf("%s: invalid phase option\n",progname);
      phase = 1;
      break;
    case 'r':
      if(sscanf(optarg,"%f",&cut.upper_theta) != 1)
	printf("%s: invalid theta option\n",progname);
      break;
    case 's':
      if(sscanf(optarg,"%d",&selector) != 1)
	printf("%s: invalid selection option\n",progname);
      if(selector == 1)
	selector = ALPHA;
      else if(selector == 2)
	selector = APERTURE;
      else
	selector = ALPHA;
      break;
    case 't':
      if(sscanf(optarg,"%d,%d,%d",&(cut.trigger[0]),&(cut.trigger[1]),
		&(cut.trigger[2])) != 3)
	printf("%s: invalid trigger option\n",progname);
      break;
    case 'u':
      if(sscanf(optarg,"%f",&cut.length_size) != 1)
	printf("%s: invalid length/size option\n",progname);
      break;
    case 'w':
      if(sscanf(optarg,"%f,%f",&cut.lower_width,&cut.upper_width) == 2){
	if(cut.upper_width < cut.lower_width){
	  tmp = cut.upper_width;
	  cut.upper_width = cut.lower_width;
	  cut.lower_width = tmp;
	}
      }else printf("%s: invalid width range option\n",progname);
      break;
    case 'y':
      if(sscanf(optarg,"%f",&(cut.lower_asymm)) != 1)
	printf("%s: invalid asymmetry option\n",progname);
      break;
    case 'z':
      if(sscanf(optarg,"%d,%d",&cut.lower_size,&cut.upper_size) != 2)
	printf("%s: invalid size option\n",progname);
      break;
    case -1:
      return;
    default:
      printf("%s: ignoring unknown option %c\n",progname,c);
      break;
    }
  }
}

/*************************************************/
/* usage                                         */
/*************************************************/  
void usage(void)
{
  printf("Usage: %s [options] gxyyyyyyp.hdf\n",progname);
  printf("       gxyyyyyyp.hdf -> parameterized file\n");
  printf("       [options]\n");
  printf("       %%d is an integer %%f is a float (real)\n");
  printf("       -a %%f      -> Alpha cut [%.2f]\n",DEFAULT_UPPER_ALPHA);
  printf("       -b %%f,%%f,%%f-> Bin size alpha,theta,image [%.2f,%.2f,%.2f]\n",
	 DEFAULT_ALPHA_BIN,DEFAULT_THETA_BIN,DEFAULT_TWOD_BIN);
  printf("       -c %%f,%%f   -> aperture Center [%.2f,%.2f]\n",
	 DEFAULT_IMAGE_XC, DEFAULT_IMAGE_YC);
  printf("       -d %%f,%%f   -> Distance cuts [%.2f,%.2f]\n",
	 DEFAULT_LOWER_DISTANCE, DEFAULT_UPPER_DISTANCE);
  printf("       -e         -> Extended supercuts [supercuts]\n");
  printf("       -f %%f      -> Frac3 cut [%.3f]\n",DEFAULT_FRAC);
  printf("       -g         -> Gaussian smoothing [TRASH CAN]\n");
  printf("       -h         -> Help\n");
  printf("       -k %%f      -> run duration [%.2lf]\n",DEFAULT_DURATION);
  printf("       -l %%f,%%f   -> Length cuts [%.2f,%.2f]\n",
	 DEFAULT_LOWER_LENGTH, DEFAULT_UPPER_LENGTH);
  printf("       -m %%f      -> sMall cuts\n",
	 DEFAULT_LENGTH_SIZE);
  printf("       -n %%d      -> Number of phase bins [%2d]\n",
	 DEFAULT_NPHASE_BINS);
  printf("       -o %%f,%%f,%%f\n");
  printf("          ,%%f,%%f,%%f,\n");
  printf("          %%f,%%f   -> Orbital binary analysis \n");
  printf("       -p %%f,%%f,%%f,\n");
  printf("          %%f,%%f,\n");
  printf("          %%d      -> Periodic analysis \n");
  printf("                     ra[radians],dec[radians],freq[Hz],\n");
  printf("                     freq'[Hz/s],epoch[MJD],\n");
  printf("                     clock[1=gps,2=osc]\n");
  printf("       -r %%f      -> apeRture size [%.2f]\n",DEFAULT_UPPER_THETA);
  printf("       -s %%d      -> Select photon data [1=alpha,2=aperture]\n");
  printf("       -t %%d,%%d,%%d-> software Trigger [%d,%d,%d]\n",
	 DEFAULT_TRIGGER_0,DEFAULT_TRIGGER_1,DEFAULT_TRIGGER_2);
  printf("       -u %%f      -> mUon cut, length/size [%.2g]\n",
	 DEFAULT_LENGTH_SIZE);
  printf("       -w %%f,%%f   -> Width cuts [%.2f,%.2f]\n",DEFAULT_LOWER_WIDTH,
	 DEFAULT_UPPER_WIDTH);
  printf("       -y %%f      -> asYmmetry cut [%.2f]\n",DEFAULT_LOWER_ASYMM);
  printf("       -z %%d,%%d   -> siZe cut [%d,%d]\n",DEFAULT_LOWER_SIZE,
	 DEFAULT_UPPER_SIZE);
  return;
}

/* SIGXY = 0.12 */
#define SIGXY2 0.0144

/*************************************************/
/* smoothing function                            */
/* 2d gaussian                                   */
/*************************************************/  
double mfun2g(double x, double y)
{
  double dx, dy;

  dx = x - (par.xo-cut.image_xc);
  dy = y - (par.yo-cut.image_yc);

  return 1.0/sqrt(2.0*M_PI*SIGXY2)*exp(-(dx*dx+dy*dy)/(2.0*SIGXY2));
}

/*************************************************/
/* smoothing function                            */
/* trash can                                     */
/*************************************************/  
double mfun2b(double x, double y)
{
  double dx, dy;

  dx = x - (par.xo-cut.image_xc);
  dy = y - (par.yo-cut.image_yc);

  if(sqrt(dx*dx + dy*dy) <= cut.upper_theta)
    return 1.0;
  else
    return 0.0;
}

/*************************************************/
/* periodic analysis                             */
/* from R. Srinivasan                            */
/*************************************************/  
void get_phase(double time, MHist *phaseh)
{
  static double deltat;
  static int iy, im, id, ny, nd, j;   
  double ph, turns, btime, otime;
  double seconds_past_midnight, t;
  static int firstpass = 0;
  
/*************************************************/
/* initial calculations                          */
/*************************************************/  
  if (firstpass == 0){
    printf("   **info** FREQ %.10lf, FREQ' %.10g at epoch %.9lf \n",
	   pulsar.freq, pulsar.freqdot, pulsar.epoch);
    printf("            RA: %.10lf DEC: %.10lf \n",parinfo.source_ra,
	   parinfo.source_dec);
    /* days since epoch (in seconds) */
    deltat = ((double)((int)(parinfo.utc_start)) - pulsar.epoch)*86400.0L; 
    iy = (int)(parinfo.idate/10000);
    im = (int)((parinfo.idate - iy*10000)/100);
    id = (int)(parinfo.idate - (iy*10000 + im*100));
    /* day of year */
    sla_calyd_(&iy,&im,&id,&ny,&nd,&j);
    if(iy > 49)
      iy += 1900;
    else
      iy += 2000;
    firstpass = 1; 
  }

/*************************************************/
/* bary needs time in seconds past midnight      */
/* from J.P. Finley                              */
/*************************************************/  
  seconds_past_midnight = (time - (double)((int)(time)))*86400.0L;
  btime = bary(seconds_past_midnight, nd, iy, parinfo.source_ra*12.0L/M_PI, 
		  parinfo.source_dec*180.0L/M_PI);
  /* solar system barycenter correction */
  par.solar_bary = btime - seconds_past_midnight;

/*************************************************/
/* binary orbital corrections                    */
/*************************************************/  
  if(orbit){
    /* get orbit uses time in MJD and returns correction in seconds */
    otime = get_orbit(time);
    par.orbit_bary = otime;
/*************************************************/
/* corrected time                                */
/*************************************************/  
    t = btime + deltat + otime;
  }else{
    t = btime + deltat;
  }

/*************************************************/
/* phase                                         */
/*************************************************/  
  ph = t*(pulsar.freq + pulsar.freqdot*0.5L*t);  
  ph = modf(ph, &turns);  
  if (ph < 0.0L) ph += 1.0L;  
  par.phase = ph;
  mfill(phaseh,ph,ph); 
  mfill(phaseh,ph+1.0L,ph+1.0L);  

  return; 
}

/*************************************************/
/* binary orbital corrections                    */
/* from T.A. Hall                                */
/*************************************************/  
double get_orbit(double time_mjd)
{
  double t, tpd;
  double orbital_phase, E, sE, cE, q, qt;
  double denom, dep;
  double omega, somega, comega, alpha, beta;

  t = time_mjd - binary.to;
  tpd = (t*86400.0L)/binary.Po;
  orbital_phase = 2.0L*M_PI*fmod(tpd - 0.5L*binary.Po_dot*tpd*tpd, 1.0L); 
  E = orbital_phase + 
    binary.ecc*sin(orbital_phase)*(1.0L + binary.ecc*cos(orbital_phase));
  denom = 1.0L - binary.ecc*cos(E);
  dep = 1.0L;

  while(fabs(dep) > 1e-12L){
    dep = (orbital_phase - (E - binary.ecc*sin(E)))/denom;
    E = E + dep;
    denom = 1.0L - binary.ecc*cos(E);
  }

  omega = binary.Omega + binary.Omega_dot*t;
  somega = sin(omega); comega = cos(omega);

  alpha = binary.asini*somega;
  beta = binary.asini*comega*sqrt(1.0L - binary.ecc*binary.ecc);
  sE = sin(E); cE = cos(E);

  q = alpha*(cE - binary.ecc) + (beta + binary.gamma)*sE;
  return(-1.0L*q + (2.0L*M_PI/binary.Po)*q*
	 (beta*cE - alpha*sE)/(1.0L - binary.ecc*cE));
}

#ifdef HAVE_MATLAB
#include "mat.h"

int MATLAB_output(void)
{
  char matname[FILENAME_MAX];
  double *value;
  int i;

  mxArray *mcut, *minfo, *mpassed;
  mxArray *field_value;
  mxArray *mtime, *msolar_bary, *morbit_bary, *menergy;
  mxArray *mphase, *mphased;
  mxArray *malpha, *mtheta, *mraw_image, *msmoothed_image;
  mxArray *msource;

  const char *cut_fieldnames[] = {"lower_size","upper_size",
				  "frac","length_size",
				  "trigger1","trigger2","trigger3",
				  "lower_width","upper_width",
				  "lower_length","upper_length",
				  "lower_distance","upper_distance",
				  "alfa","asymm",
				  "theta","center_x","center_y"};
  const char *info_fieldnames[] = {"idate","utc_start","source_ra",
				   "source_dec","duration"};

  MATFile *matfile;

  // JPF - I am putting in this change as per SJF in ver. 3.1.0 to fix
  // the MatLab problems
  // SJF - this next macro copies the memory beloning to a MATLAB matrix
  // into another MATLAB matrix. This must be done since MATLAB wants to
  // free both of them when they are finished -- you cannot just assign
  // the pointer as RWL had been doing!

  // 1-d matrix (vector)
#define COPYMBOOK1(A,M) memcpy(mxGetPr(A),(M).matrixptr,(M).nibins*sizeof(double));

  // 2-d matrix
#define COPYMBOOK2(A,M) memcpy(mxGetPr(A),(M).matrixptr,(M).nibins*(M).njbins*sizeof(double));

#define COPYMBOOK3(A,M,N)   memcpy(mxGetPr(A), M, N*sizeof(double));

  strncpy(matname,(char *)basename(filename),8); matname[8] = '\0';
  strncat(matname,".mat",4);

  minfo = mxCreateStructMatrix(1,1,5,info_fieldnames);
  mxSetName(minfo,"info");

  field_value = mxCreateDoubleMatrix(1,1,mxREAL);
  *mxGetPr(field_value) = (double)parinfo.idate;
  mxSetField(minfo,0,"idate",field_value);

  field_value = mxCreateDoubleMatrix(1,1,mxREAL);
  *mxGetPr(field_value) = parinfo.utc_start;
  mxSetField(minfo,0,"utc_start",field_value);

  field_value = mxCreateDoubleMatrix(1,1,mxREAL);
  *mxGetPr(field_value) = parinfo.source_ra;
  mxSetField(minfo,0,"source_ra",field_value);

  field_value = mxCreateDoubleMatrix(1,1,mxREAL);
  *mxGetPr(field_value) = parinfo.source_dec;
  mxSetField(minfo,0,"source_dec",field_value);

  field_value = mxCreateDoubleMatrix(1,1,mxREAL);
  *mxGetPr(field_value) = parinfo.duration;
  mxSetField(minfo,0,"duration",field_value);

  mcut = mxCreateStructMatrix(1,1,18,cut_fieldnames);
  mxSetName(mcut,"cut");

  field_value = mxCreateDoubleMatrix(1,1,mxREAL);
  *mxGetPr(field_value) = (double)cut.lower_size; 
  mxSetField(mcut,0,"lower_size",field_value);

  field_value = mxCreateDoubleMatrix(1,1,mxREAL);
  *mxGetPr(field_value) = (double)cut.upper_size;
  mxSetField(mcut,0,"upper_size",field_value);

  field_value = mxCreateDoubleMatrix(1,1,mxREAL);
  *mxGetPr(field_value) = (double)cut.frac; 
  mxSetField(mcut,0,"frac",field_value);

  field_value = mxCreateDoubleMatrix(1,1,mxREAL);
  *mxGetPr(field_value) = (double)cut.length_size;
  mxSetField(mcut,0,"length_size",field_value);


  field_value = mxCreateDoubleMatrix(1,1,mxREAL);
  *mxGetPr(field_value) = (double)cut.trigger[0]; 
  mxSetField(mcut,0,"trigger1",field_value);


  field_value = mxCreateDoubleMatrix(1,1,mxREAL);
  *mxGetPr(field_value) = (double)cut.trigger[1];
  mxSetField(mcut,0,"trigger2",field_value);


  field_value = mxCreateDoubleMatrix(1,1,mxREAL);
  *mxGetPr(field_value) = (double)cut.trigger[2]; 
  mxSetField(mcut,0,"trigger3",field_value);

  field_value = mxCreateDoubleMatrix(1,1,mxREAL);
  *mxGetPr(field_value) = (double)cut.lower_width; 
  mxSetField(mcut,0,"lower_width",field_value);

  field_value = mxCreateDoubleMatrix(1,1,mxREAL);
  *mxGetPr(field_value) = (double)cut.upper_width;
  mxSetField(mcut,0,"upper_width",field_value);

  field_value = mxCreateDoubleMatrix(1,1,mxREAL);
  *mxGetPr(field_value) = (double)cut.lower_length; 
  mxSetField(mcut,0,"lower_length",field_value);

  field_value = mxCreateDoubleMatrix(1,1,mxREAL);
  *mxGetPr(field_value) = (double)cut.upper_length;
  mxSetField(mcut,0,"upper_length",field_value);

  field_value = mxCreateDoubleMatrix(1,1,mxREAL);
  *mxGetPr(field_value) = (double)cut.lower_distance;
  mxSetField(mcut,0,"lower_distance",field_value);

  field_value = mxCreateDoubleMatrix(1,1,mxREAL);
  *mxGetPr(field_value) = (double)cut.upper_distance;
  mxSetField(mcut,0,"upper_distance",field_value);

  field_value = mxCreateDoubleMatrix(1,1,mxREAL);
  *mxGetPr(field_value) = (double)cut.upper_alpha;
  mxSetField(mcut,0,"alfa",field_value);

  field_value = mxCreateDoubleMatrix(1,1,mxREAL);
  *mxGetPr(field_value) = (double)cut.lower_asymm;
  mxSetField(mcut,0,"asymm",field_value);

  field_value = mxCreateDoubleMatrix(1,1,mxREAL);
  *mxGetPr(field_value) = (double)cut.upper_theta;
  mxSetField(mcut,0,"theta",field_value);

  field_value = mxCreateDoubleMatrix(1,1,mxREAL);
  *mxGetPr(field_value) = (double)cut.image_xc;
  mxSetField(mcut,0,"center_x",field_value);

  field_value = mxCreateDoubleMatrix(1,1,mxREAL);
  *mxGetPr(field_value) = (double)cut.image_yc;
  mxSetField(mcut,0,"center_y",field_value);

  mpassed = mxCreateDoubleMatrix(5,1,mxREAL);
  mxSetName(mpassed,"passed");
  value = mxGetPr(mpassed);
  for(i=0;i < 5;i++)
    value[i] = (double)output_results.passed[i]; 
    
  msource = mxCreateString(parinfo.source_name);
  mxSetName(msource,"source");

  mtime = mxCreateDoubleMatrix(output_results.passed[selector],1,mxREAL);
  mxSetName(mtime,"time");
  // JPF&MC This is necessary for ver. 4 since the structures have changed
  // memcpy(mxGetPr(mtime), output_results.photon.time, 100000*sizeof(double));
  COPYMBOOK3(mtime,output_results.photon.time,output_results.passed[selector]);
  // JPF WRONG mxSetPr(mtime,output_results.photon.time);

  msolar_bary = mxCreateDoubleMatrix(output_results.passed[selector],1,mxREAL);
  mxSetName(msolar_bary,"solar_bary"); 
  COPYMBOOK3(msolar_bary,output_results.photon.solar_bary,output_results.passed[selector]);
  // COPYMBOOK1(msolar_bary,output_results.photon.solar_bary);
  // JPF WRONG mxSetPr(msolar_bary,output_results.photon.solar_bary);
    
  morbit_bary = mxCreateDoubleMatrix(output_results.passed[selector],1,mxREAL);
  mxSetName(morbit_bary,"orbit_bary");
  COPYMBOOK3(morbit_bary,output_results.photon.orbit_bary,output_results.passed[selector]);
  // COPYMBOOK1(morbit_bary,output_results.photon.orbit_bary);
  // JPF WRONG mxSetPr(morbit_bary,output_results.photon.orbit_bary);

  mphase = mxCreateDoubleMatrix(output_results.passed[selector],1,mxREAL);
  mxSetName(mphase,"phase");
  COPYMBOOK3(mphase,output_results.photon.phase,output_results.passed[selector]);
  // COPYMBOOK1(mphase,output_results.photon.phase);
  // JPF WRONG mxSetPr(mphase,output_results.photon.phase);

  menergy = mxCreateDoubleMatrix(output_results.passed[selector],1,mxREAL);
  mxSetName(menergy,"energy");
  COPYMBOOK3(menergy,output_results.photon.energy,output_results.passed[selector]);
  // COPYMBOOK1(menergy,output_results.photon.energy);
  // JPF WRONG mxSetPr(menergy,output_results.photon.energy);

  mphased = mxCreateDoubleMatrix(output_results.phase.nibins,1,mxREAL);
  mxSetName(mphased,"phased"); 
  COPYMBOOK1(mphased,output_results.phase);
  // JPF WRONG mxSetPr(mphased,output_results.phase.matrixptr);

  malpha = mxCreateDoubleMatrix(output_results.alpha.nibins,1,mxREAL);
  mxSetName(malpha,"alfa"); 
  COPYMBOOK1(malpha,output_results.alpha);
  // JPF WRONG mxSetPr(malpha,output_results.alpha.matrixptr);

  mtheta = mxCreateDoubleMatrix(output_results.theta.nibins,1,mxREAL);
  mxSetName(mtheta,"theta"); 
  COPYMBOOK1(mtheta,output_results.theta);
  // JPF WRONG mxSetPr(mtheta,output_results.theta.matrixptr);

  mraw_image = mxCreateDoubleMatrix(output_results.raw_image.nibins,
				output_results.raw_image.njbins,mxREAL);
  mxSetName(mraw_image,"raw_image"); 
  COPYMBOOK2(mraw_image,output_results.raw_image);
  // JPF WRONG mxSetPr(mraw_image,output_results.raw_image.matrixptr);

  msmoothed_image = mxCreateDoubleMatrix(output_results.smoothed_image.nibins,
				 output_results.smoothed_image.njbins,mxREAL);
  mxSetName(msmoothed_image,"smoothed_image"); 
  COPYMBOOK2(msmoothed_image,output_results.smoothed_image);
  // JPF WRONG mxSetPr(msmoothed_image,output_results.smoothed_image.matrixptr);

  matfile = matOpen(matname,"w");

  matPutArray(matfile,mcut);
  matPutArray(matfile,minfo);
  matPutArray(matfile,msource);
  matPutArray(matfile,mpassed);  

  matPutArray(matfile,mtime);
  matPutArray(matfile,msolar_bary); 
  matPutArray(matfile,morbit_bary);
  matPutArray(matfile,mphase);
  matPutArray(matfile,menergy);
  
  matPutArray(matfile,mphased);
  matPutArray(matfile,malpha); 
  matPutArray(matfile,mtheta);
  matPutArray(matfile,mraw_image); 
  matPutArray(matfile,msmoothed_image);

  matClose(matfile);

  mxDestroyArray(minfo);
  mxDestroyArray(mcut);
  mxDestroyArray(msource);
  mxDestroyArray(mpassed); 
  
  mxDestroyArray(mphased);
  mxDestroyArray(malpha); 
  mxDestroyArray(mtheta); 
  mxDestroyArray(mraw_image); 
  mxDestroyArray(msmoothed_image);

  return 1;
}
#endif
