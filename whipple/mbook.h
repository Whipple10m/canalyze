/*************************************************/
/* storage for histograms                        */
/*************************************************/
struct matrix_s{
  double *matrixptr;
  double imin;
  double jmin;
  double imax;
  double jmax;
  double istep;
  double jstep;
  double *imid;
  double *jmid;
  int nibins;
  int njbins;
};
typedef struct matrix_s MHist;

/*************************************************/
/* Function declarations                         */
/*************************************************/
void mbook1(MHist *, int , double, double);
void mbook2(MHist *, int, double, double, int, double, double);
int mfill(MHist *, double, double);
int mfun(MHist *, double (*)(double, double));
