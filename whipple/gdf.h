/**************************************************************************/
/* definition of all records available to main program                    */
/**************************************************************************/
#ifdef SUN
#define short int
#endif

#define gdf_version 83
/**************************************************************************/
/*                  global parameter                                      */
/**************************************************************************/

#define gdf_tele_10  1   /* 10 m            */
#define gdf_tele_11  2   /* 11 m            */
#define gdf_tele_max 2   /* # of telescopes */

/**************************************************************************/
/*                  sky quality                                           */
/**************************************************************************/

#define sky_good      1   /* definitely good data (A) */
#define sky_uncertain 2   /* probably good data (B)   */
#define sky_bad       3   /* definitely not good (C)  */

/**************************************************************************/
/*                  trigger mode                                          */
/**************************************************************************/

#define trig_mode_single 1   /* single threshold */
#define trig_mode_dual   2   /* double threshold */

/* The old trigger bit position are arranged according to the GRALP
   trigger numbers. When old GRALP files are converted to this format
   the trigger bit corresponding to the event code is set. */


#define gdf_trig_short 1   /* GRALP only, obsolete                      */
#define gdf_trig_long  2   /* GRALP only, obsolete                      */
#define gdf_trig_test1 3   /* fixed period trigger, no TDC time delay   */
#define gdf_trig_test2 4   /* fixed period trigger, with TDC time delay */
#define gdf_trig_wwvb  6   /* WWVB time marker                          */
#define gdf_trig_stime 7   /* siderial time marker, obsolete            */
#define gdf_trig_hig   8   /* high level trigger                        */
#define gdf_trig_low   9   /* low level trigger                         */
#define gdf_trig_eas   12  /* GRALP, EAS trigger, obsolete              */

/* The Hytec based system has its own trigger bits conventions.
   Note that a given event may have more than one trigger bit set. */

#define gdf_trig_ped   0   /* pedestal trigger                          */
#define gdf_trig_pst   1   /* pattern selection trigger                 */
#define gdf_trig_mul   2   /* multiplicity trigger                      */

/**************************************************************************/
/*                  type of run                                           */
/**************************************************************************/

#define gdf_run_type_stereo 1
#define gdf_run_type_10     2
#define gdf_run_type_11     3
#define gdf_run_type_mc     4

/**************************************************************************/
/*                  telescope tracking status/mode                        */
/**************************************************************************/

#define gdf_track_mode_on      1
#define gdf_track_mode_off     2
#define gdf_track_mode_slewing 3
#define gdf_track_mode_standby 4
#define gdf_track_mode_zenith  5
#define gdf_track_mode_check   6   /* pointing check    */
#define gdf_track_mode_stowing 7   /* stowing telescope */
#define gdf_track_mode_drift   8   /* drift scan        */
#define gdf_track_mode_max     8

/**************************************************************************/
/*                  run information                                       */
/**************************************************************************/

#define gdf_run_mcl 16000   /* max. comment length in characters */

struct gdf_run_s {
  int      version;                     /* software version number           */
  int      reserved;                    /* for future use                    */
  int      checksum[2];                 /* data part check sum               */
  double   utc;                         /* current UTC [mjd]                 */
  unsigned status[gdf_tele_max];        /* detector status [bits]            */
  int      idate;                       /* VAX date, local time [yymmdd]     */
  int      itime;                       /* VAX time, local time [hhmmss]     */
  int      year;                        /* Gregorian year, UTC time          */
  int      run;                         /* run number                        */
  int      type;                        /* type of run                       */
  int      sky_quality;                 /* sky quality                       */
  int      trig_mode[gdf_tele_max];     /* trigger setup                     */
  int      trig_nhig[gdf_tele_max];     /* min tubes above high threshold    */
  int      trig_nlow[gdf_tele_max];     /* min tubes above low threshold     */
  int      clen;                        /* actual length of comment [bytes]  */
  float    sid_length[1];               /* siderial nominal run length [min] */
  float    sid_cycle;                   /* siderial nominal cycle time [min] */
  float    sid_actual;                  /* siderial actual logged time [min] */
  float    trig_thrlow[gdf_tele_max];   /* low level trigger thresholds [V]  */
  float    trig_thrhig[gdf_tele_max];   /* high level trigger thresholds [V] */
  double   utc_start;                   /* nominal UTC start of run [mjd]    */
  double   utc_end;                     /* nominal UTC end of run [mjd]      */
  char     file[80];                    /* file name                         */
  char     observer[80];                /* observer names                    */
  char     comment[gdf_run_mcl];        /* any observer comments             */
  unsigned new;                         /* TRUE if just read from file       */
  unsigned valid;                       /* TRUE if record contents valid     */
};

/**************************************************************************/
/*                       CCD camera results                               */
/**************************************************************************/

#define gdf_ccd_nstar_max 100   /* max number of stars stored */

struct gdf_ccd_s {
  int version;                       /* software version number         */
  int reserved;                      /* for future use                  */
  int checksum[2];                   /* data part check sum             */
  double utc;                        /* UTC time picture finished [mjd] */
  int telescope;                     /* telescope identifier            */
  int nstar;                         /* actual number of stars          */
  int cycle;                         /* number of updates sofar         */
  int interlace;                     /* interlace: on=1,off=0           */
  int bias;                          /* pedestal                        */
  int gain;                          /* gain                            */
  int noise_range;                   /*                                 */
  int exposure;                      /* exposure time [msec]            */
  int intensity[gdf_ccd_nstar_max];  /* intensity                       */
  int status[gdf_ccd_nstar_max];     /* status                          */
  float marker[2];                   /* marker position of center pmt   */
  float dark_mean;                   /* dark file mean                  */
  float dark_sigma;                  /* dark file sigma                 */
  float low_mean;                    /* mean, low live pixel            */
  float low_sigma;                   /* sigma, low live pixel           */
  float noise_threshold;             /*                                 */
  float noise_level;                 /*                                 */
  float star[2*gdf_ccd_nstar_max];   /* 1=x,2=y                         */
  unsigned new;                      /* TRUE if just read from file     */
  unsigned valid;                    /* TRUE if record contents valid   */
};

/**************************************************************************/
/*                       Tracking                                         */
/**************************************************************************/

struct gdf_track_s {
  int version;          /* software version number                    */
  int reserved;         /* for future use                             */
  int checksum[2];      /* data part check sum                        */
  double utc;           /* current UTC time [mjd]                     */
  int telescope;        /* telescope id                               */
  int mode;             /* tracking mode                              */
  int cycle;            /* number of updates sofar                    */
  unsigned status;      /* telescope status [bits]                    */
  double rasc_2000;     /* source right ascension, FK5 J2000 [rad]    */
  double decl_2000;     /* source declination, FK5 J2000 [rad]        */
  double rasc_today;    /* source right ascension, FK5 today [rad]    */
  double decl_today;    /* source declination, FK5 today [rad]        */
  double rasc_tele;     /* telescope right ascension, FK5 today [rad] */
  double decl_tele;     /* telescope declination, FK5 today [rad]     */
  double azimuth;       /* telescope pointing, +west, north = 0 [rad] */
  double elevation;     /* telescope pointing, +up, horizon = 0 [rad] */
  double deviation;     /* angle nominal/actual position [rad]        */
  double rasc_offset;   /* RA offset for off-source runs              */
  double decl_offset;   /* DE offset for off-source runs              */
  double stl;           /* local siderial time [rad]                  */
  double height;        /* height of interaction point                */
  double azi_incl;      /* azimuth change for inclination [rad]       */
  double ele_incl;      /* elevation change for incl. [rad]           */
  char source[80];      /* source name                                */
  unsigned new;         /* TRUE if just read from file                */
  unsigned valid;       /* TRUE if record contents valid              */
};

/**************************************************************************/
/*                  high voltage                                          */
/**************************************************************************/

#define gdf_hvc_max 640   /* max number of HV channels, (16 per module, 11m) */

struct gdf_hv_s {
  int version;                 /* software version number               */
  int reserved;                /* for future use                        */
  int checksum[2];             /* data part check sum                   */
  double utc;                  /* UTC time, end of measurement [hhmmss] */
  int telescope;               /* telescope identifier                  */
  int mode;                    /* current operation mode                */
  int nch;                     /* channels voltages values              */
  int cycle;                   /* read cycle number                     */
  short status[gdf_hvc_max];   /* status of each HV channel [bits]      */
  float v_set[gdf_hvc_max];    /* presently set voltage [V]             */
  float v_actual[gdf_hvc_max]; /* actual measured voltage [V]           */
  float i_supply[gdf_hvc_max]; /* HV supply current [uA]                */
  float i_anode[gdf_hvc_max];  /* measured anode current [uA]           */
  unsigned new;                /* TRUE if just read from file           */
  unsigned valid;              /* TRUE if record contents valid         */
};

/**************************************************************************/
/*                  detector calibration: peds and gains                  */
/**************************************************************************/

#define gdf_cal_npm_max 541   /* max number of tubes */

struct gdf_cal_s {
  int version;                             /* software version number       */
  int reserved;                            /* for future use                */
  int checksum[2];                         /* data part of check sum        */
  double utc;                              /* UTC time first event [mjd]    */
  int npm;                                 /* number of photomultipliers    */
  int run;                                 /* run number                    */
  int method;                              /* method used in calculation    */
  int status[gdf_cal_npm_max];             /* tube status                   */
  float peak[gdf_cal_npm_max];             /* pedestal maximum              */
  float width[gdf_cal_npm_max];            /* pedestal width                */
  float pedestal[gdf_cal_npm_max];         /* value itself                  */
  float asycor[gdf_cal_npm_max];           /* ped linear width correction   */
  float symcor[gdf_cal_npm_max];           /* ped quadr. width correction   */
  float gain[gdf_cal_npm_max];             /* gain correction factor        */
  float exponent[gdf_cal_npm_max];         /* power law gain corr.          */
  float peak_error[gdf_cal_npm_max];       /* error of peak                 */
  float width_error[gdf_cal_npm_max];      /* error of width                */
  float pedestal_error[gdf_cal_npm_max];   /* error of pedestal             */
  float asycor_error[gdf_cal_npm_max];     /* error of linear correction    */
  float symcor_error[gdf_cal_npm_max];     /* error of quadr. correction    */
  float gain_error[gdf_cal_npm_max];       /* error of gain correction      */
  float exponent_error[gdf_cal_npm_max];   /* error of power law corr.      */
  float pedestal_pro[gdf_cal_npm_max];     /* probability pedestal ok       */
  float gain_pro[gdf_cal_npm_max];         /* probability gain ok           */
  unsigned new;                            /* TRUE if just read from file   */
  unsigned valid;                          /* TRUE if record contents valid */
};

/**************************************************************************/
/*                  10 meter frame                                        */
/**************************************************************************/

#define gdf_fr10_nadc 636   /* max # of ADCs (12 per module)       */
#define gdf_fr10_nsca 640   /* max # of phase TDCs (32 per module) */ 
#define gdf_fr10_nphs 8     /* max # of scaler (8 per module)      */

struct gdf_fr10_s {
  int version;                     /* software version number              */
  int reserved;                    /* for future use                       */
  int checksum[2];                 /* data part check sum                  */
  double utc;                      /* current UTC time [mjd]               */
  unsigned status;                 /* detector status [bits]               */
  unsigned mark_gps;               /* last recorded GPS (perhaps!) [50ns]  */
  int nphs;                        /* number of phase TDCs                 */
  int nadc;                        /* number of ADCs                       */
  int nsca;                        /* number of scalers                    */
  int run;                         /* run number                           */
  int frame;                       /* frame number                         */
  int gps_mjd;                     /* soon: modified Julian days [mjd]     */
  int gps_sec;                     /* soon: seconds since midnight [sec]   */
  int gps_ns;                      /* soon: time from last GPS second [ns] */
  short cal_adc[gdf_fr10_nadc];    /* ADC with internal test voltage       */
  short ped_adc1[gdf_fr10_nadc];   /* ADC first random event               */
  short ped_adc2[gdf_fr10_nadc];   /* ADC second random event              */
  short scalc[gdf_fr10_nsca];      /* current monitor scaler               */
  short scals[gdf_fr10_nsca];      /* single rates scaler                  */
  short gps_clock[3];              /* GPS time last event [bits]           */
  short phase_delay;               /* phase delay module settings          */
  short phs1[gdf_fr10_nphs];       /* phase TDC first random event         */
  short phs2[gdf_fr10_nphs];       /* phase TDC second random event        */
  short gps_status[2];             /* soon: GPS status flags               */
  int align;                       /* # 64 bit alignment                   */
  unsigned new;                    /* TRUE if just read from file          */
  unsigned valid;                  /* TRUE if record contents valid        */
};

/**************************************************************************/
/*                  11 meter frame                                        */
/**************************************************************************/

#define gdf_fr11_nadc 180   /* max number of ADCs (12 per module)      */
#define gdf_fr11_ntdc 176   /* max number of TDCs (8 per module)       */
#define gdf_fr11_nphs 8     /* max number of phase TDCs (8 per module) */
#define gdf_fr11_nsca 192   /* max number of scalers (32 per module)   */

struct gdf_fr11_s {
  int version;                     /* software version number              */
  int reserved;                    /* for future use                       */
  int checksum[2];                 /* data part check sum                  */
  double utc;                      /* current UTC time [mjd]               */
  unsigned status;                 /* detector status [bits]               */
  unsigned mark_gps;               /* last recorded GPS (perhaps!) [50ns]  */
  int nphs;                        /* number of phase TDCs                 */
  int ntdc;                        /* number of TDCs                       */
  int nadc;                        /* number of ADCs                       */
  int nsca;                        /* number of scalers                    */
  int run;                         /* run number                           */
  int frame;                       /* frame number                         */
  int gps_mjd;                     /* soon: modified Julian days [mjd]     */
  int gps_sec;                     /* soon: seconds since midnight [sec]   */
  int gps_ns;                      /* soon: time from last GPS second [ns] */
  int align_1;                     /* # 64 bit alignment                   */
  short cal_adc[gdf_fr11_nadc];    /* ADC with internal test voltage       */
  short ped_adc1[gdf_fr11_nadc];   /* ADC first random event               */
  short ped_adc2[gdf_fr11_nadc];   /* ADC second random event              */
  short tdc1[gdf_fr11_ntdc];       /* ADC first random event               */
  short tdc2[gdf_fr11_ntdc];       /* ADC second random event              */
  short scalc[gdf_fr11_nsca];      /* current monitor scaler               */
  short scals[gdf_fr11_nsca];      /* single rates scaler                  */
  short geos_clock[3];             /* Geos time last event [bits]          */
  short phase_delay;               /* phase delay module settings          */
  short gps_status[2];             /* soon: GPS status flags               */
  short phs1[gdf_fr11_nphs];       /* phase TDC first random event         */
  short phs2[gdf_fr11_nphs];       /* phase TDC second random event        */
  int align_2;                     /* # 64 bit alignment                   */
  unsigned new;                    /* TRUE if just read from file          */
  unsigned valid;                  /* TRUE if record contents valid        */
};

/**************************************************************************/
/*                  10 meter event                                        */
/**************************************************************************/

#define gdf_ev10_nadc  636   /* max # of ADCs (12 per module)      */
#define gdf_ev10_nphs  8     /* max # of phase TDCs (8 per module) */
#define gdf_ev10_nbrst 12    /* max number of burst scalers        */
#define gdf_ev10_ntrg  65    /* pattern trigger data words         */

struct gdf_ev10_s {
  int version;                   /* software version number              */
  int reserved;                  /* # for future use                     */
  int checksum[2];               /* # data part check sum                */
  double utc;                    /* GPS UTC time of event [mjd]          */
  int nadc;                      /* number of ADCs                       */
  int run;                       /* run number                           */
  int event;                     /* event number                         */
  int live_sec;                  /* live time from start of run [sec]    */
  int live_ns;                   /* live time last incomplete sec [ns]   */
  int frame;                     /* frame number                         */
  int frame_event;               /* events within frame                  */
  int abort_cnt;                 /* number of aborts in frame            */
  int nphs;                      /* number of phase TDCs                 */
  int nbrst;                     /* number of burst scalers              */
  int gps_mjd;                   /* soon: modified Julian days [mjd]     */
  int gps_sec;                   /* soon: seconds since midnight [sec]   */
  int gps_ns;                    /* soon: time from last GPS second [ns] */
  int ntrg;                      /* number of trigger patterns           */
  int elapsed_sec;               /* sec from start of run [sec]          */
  int elapsed_ns;                /* ns since last elapsed sec            */
  int grs_clock[3];              /* Wisconsin TrueTime interface         */
  int align;                     /* 64 bit alignment                     */
  unsigned trigger;              /* trigger bit pattern                  */
  unsigned status;               /* detector status flags [bits]         */
  unsigned mark_gps;             /* last recorded GPS [50ns]             */
  unsigned mark_open;            /* last calibration mark open [50ns]    */
  unsigned mark_close;           /* last calibration mark close [50ns]   */
  unsigned gate_open;            /* last event gate open [50ns]          */
  unsigned gate_close;           /* last event gate close [50ns]         */
  unsigned pattern[gdf_ev10_ntrg];/*pattern trigger output               */
  short adc[gdf_ev10_nadc];      /* event ADC's                          */
  short gps_clock[3];            /* GPS satellite time [bits]            */
  short phase_delay;             /* phase delay module settings          */
  short phs[gdf_ev10_nphs];      /* phase TDCs                           */
  short burst[gdf_ev10_nbrst];   /* burst scalers                        */
  short gps_status[2];           /* soon: GPS status flags               */
  short track[2];                /* telscope angle encoders [bits]       */
  unsigned new;                  /* TRUE if just read from file          */
  unsigned valid;                /* TRUE if record contents valid        */
};

/**************************************************************************/
/*                  11 meter event                                        */
/**************************************************************************/

#define gdf_ev11_nadc  180   /* ADCs (12 per module) */
#define gdf_ev11_ntdc  176   /* TDCs (8 per module)  */
#define gdf_ev11_nphs  8     /* phase TDCs           */
#define gdf_ev11_nbrst 12    /* burst scalers        */

struct gdf_ev11_s {
  int version;                   /* software version number              */
  int reserved;                  /* # for future use                     */
  int checksum[2];               /* # data part check sum                */
  double utc;                    /* GPS UTC time of event [mjd]          */
  unsigned trigger;              /* trigger bit pattern                  */
  unsigned status;               /* detector status flags [bits]         */
  unsigned mark_gps;             /* last recorded GPS [50ns]             */
  unsigned mark_open;            /* last calibration mark open [50ns]    */
  unsigned mark_close;           /* last calibration mark close [50ns]   */
  unsigned gate_open;            /* last event gate open [50ns]          */
  unsigned gate_close;           /* last event gate close [50ns]         */
  int nbrst;                     /* number of burst scalers              */
  int nphs;                      /* number of phase TDCs                 */
  int ntdc;                      /* number of TDCs                       */
  int nadc;                      /* number of ADCs                       */
  int run;                       /* run number                           */
  int event;                     /* event number                         */
  int live_sec;                  /* live time from start of run [sec]    */
  int live_ns;                   /* live time last incomplete sec [ns]   */
  int frame;                     /* frame number                         */
  int frame_event;               /* events within frame                  */
  int abort_cnt;                 /* number of aborts in frame            */
  int gps_mjd;                   /* soon: modified Julian days [mjd]     */
  int gps_sec;                   /* soon: seconds since midnight [sec]   */
  int gps_ns;                    /* soon: time from last GPS second [ns] */
  int align;                     /* # 64 bit alignment                   */
  short adc[gdf_ev11_nadc];      /* event ADC's                          */
  short tdc[gdf_ev11_ntdc];      /* event TDC's                          */
  short geos_clock[3];           /* Geos satellite time [bits]           */
  short phase_delay;             /* phase delay module settings          */
  short gps_status[2];           /* soon: GPS status flags               */
  short phs[gdf_ev11_nphs];      /* phase TDCs                           */
  short burst[gdf_ev11_nbrst];   /* burst scalers                        */
  short track[2];                /* telscope angle encoders [bits]       */
  unsigned new;                  /* TRUE if just read from file          */
  unsigned valid;                /* TRUE if record contents valid        */
};

/**************************************************************************/
/*                  monte carlo event header                              */
/**************************************************************************/

#define gdf_mirror_max 10
#define gdf_mce_r 3+16*(5+8*gdf_mirror_max)
#define gdf_mce_i 2+16*3

struct gdf_mce_s {
  int version;                        /* software version number           */
  int reserved;                       /* # for future use                  */
  int checksum[2];                    /* # data part check sum             */
  double utc;                         /* UTC time of event [mjd]           */
  int charge;                         /* incident particle: charge         */
  int nucleons;                       /* incident particle: nucleons       */
  int mirrors;                        /* total number of mirrors           */
  float momentum;                     /* incident particle: momentum       */
  float elevation;                    /* incident par:el.(up=90) [deg]     */
  float azimuth;                      /* incident par:az.(E=0,N=90) [deg]  */
  float height;                       /* height above sea level [m]        */
  float radius[gdf_mirror_max];       /* mirror radius [m]                 */
  float position[3*gdf_mirror_max];   /* geographic {x,y,z} [m]            */
  float location[3*gdf_mirror_max];   /* shower coordinates {x',y',z'} [m] */
  float time[gdf_mirror_max];         /* relative arrival time [ns]        */
  float align;                        /* dummy word for 64-bit alignment   */
  unsigned new;                       /* TRUE if just read from file       */
  unsigned valid;                     /* TRUE if record contents valid     */
};
  
/**************************************************************************/
/*     monte carlo photons hitting individual mirrors                     */
/**************************************************************************/

#define gdf_photon_max 20000   /* max. possible # of photons */

struct gdf_photon_s {
  float xu;       /* x unit velocity (shower)        */
  float yu;       /* y unit velocity (shower)        */
  float xi;       /* horiz impact point (shower) [m] */
  float yi;       /* vert impact point (shower) [m]  */
  float xf;       /* x-pos focal plane [m]           */
  float yf;       /* y-pos focal plane [m]           */
  float energy;   /* photon energy [eV]              */
  float height;   /* emission height (geo) [m]       */
  float time;     /* relative arrival time [ns]      */
  float weight;   /* Monte Carlo weight              */
};

struct gdf_mcp_s {
  int version;                                  /* software version number */
  int reserved;                                 /* # for future use        */
  int checksum[2];                              /* # data part check sum   */
  double dummy;                                 /* would normally be UTC   */
  int photons;                                  /* actual # of photons     */
  int mirror[gdf_photon_max];                   /* mirror number           */
  struct gdf_photon_s photon[gdf_photon_max];   /* use this structure for
                                                   this many photons       */
  float align;                                  /* #dummy word for 64-bit
                                                   alignment               */
  unsigned new;                                 /* TRUE if just read from
                                                   file                    */
  unsigned valid;                               /* TRUE if record contents 
                                                   valid                   */
};

/**************************************************************************/
/*   the C equivalent of the FORTRAN common block, note the underscore    */
/*   after common block name.                                             */
/**************************************************************************/
  
struct {
  int gdf_data_first;                /* first 32-bit word       */
  int gdf_data_align;                /* to get 64-bit alignment */
  struct gdf_run_s gdf_run;
  struct gdf_ccd_s gdf_ccd[2];
  struct gdf_track_s gdf_track[2];
  struct gdf_hv_s gdf_hv[2];
  struct gdf_fr10_s gdf_fr10;
  struct gdf_ev10_s gdf_ev10;
  struct gdf_fr11_s gdf_fr11;
  struct gdf_ev11_s gdf_ev11;
  struct gdf_mce_s gdf_mce;
  struct gdf_mcp_s gdf_mcp;
  struct gdf_cal_s gdf_cal[2];
  int gdf_data_last;
} gdf_data_;
