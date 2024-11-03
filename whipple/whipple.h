#define NNEIGHBOURS     6
#define RADTODEG        57.29577951

struct coords {
  float *wxdeg;
  float *wydeg;
};

/* General function declarations */
struct coords cwhipplecoords();
int ** cwhippleneighbours();
int getpeds();
int getnlist();
int getn2gains();
int getpair();
void ggetgtoff();
int get_npmt(int);
float rad_to_ddmmss(float);
float rad_to_hhmmss(float);
float hhmmss_to_rad(float);
float ddmmss_to_rad(float);

