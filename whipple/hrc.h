struct header {
  char   *runid, *mode, *sourcename, *skyq, *comms;
  int   nevents, kdur, kdate, kut, kst;
  double duration, mjd, frjd, gpsbeg;
  float  ra, dec, azimuth, elevation;
  int flag;
};
struct frame {
  int   code;
  short *channel;
  double time;
  int flag;
};
  
struct header get_head();
struct frame get_frame();
