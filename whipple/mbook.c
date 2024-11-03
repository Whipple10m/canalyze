/****************************************************************************/
/* mbook: histogram data storage routines                                   */
/* Written by R.W. Lessard 961214                                           */
/****************************************************************************/
#include <stdio.h>
#include <math.h>
#include "mbook.h"

/****************************************************************************/
/* mfill: fills histogram with data, checks bounds                          */
/****************************************************************************/
int mfill(MHist *matrix, double ivalue, double jvalue)
{
  int ibin, jbin;

  if(matrix->njbins == 0){ /* 1-dimensional histogram */
    if((ivalue <= matrix->imax)&&(ivalue >= matrix->imin)){
      ibin = (int)((ivalue - matrix->imin)/matrix->istep);
      if((ibin < matrix->nibins) && (ibin >= 0))
	*(matrix->matrixptr + ibin)+=1.0;
      return 1;
    }else{
      return 0;
    }
  }
  if((ivalue <= matrix->imax)&&(ivalue >= matrix->imin)&&
     (jvalue <= matrix->jmax)&&(jvalue >= matrix->jmin)){ /* 2-d */
    ibin = (int)((ivalue - matrix->imin)/matrix->istep);
    jbin = (int)((jvalue - matrix->jmin)/matrix->jstep);
    if((ibin < matrix->nibins) && (jbin < matrix->njbins))
      *(matrix->matrixptr + ibin + jbin*matrix->nibins)+=1.0;
    return 1;
  }
  return 0;
}

/****************************************************************************/
/* mfun: fills histogram with data from external function                   */
/****************************************************************************/
int mfun(MHist *matrix, double (*func)(double, double))
{
  int i, j;

  if(matrix->njbins == 0){
    for(i=0;i < matrix->nibins;i++)
      *(matrix->matrixptr + i)+=func(*(matrix->imid+i),*(matrix->imid+i));
  }
  for(i=0;i < matrix->nibins;i++){
    for(j=0;j < matrix->njbins;j++)
      *(matrix->matrixptr + i + j*matrix->nibins)+=
	func(*(matrix->imid+i),*(matrix->jmid+j));
  }
  return 1;
}

/****************************************************************************/
/* mbook1: set up 1-d histogram                                             */
/****************************************************************************/
void mbook1(MHist *matrix, int nbins, double min, double max)
{
  int i;

  matrix->matrixptr = (double *)malloc(nbins*sizeof(double));
  matrix->imid = (double *)malloc(nbins*sizeof(double));
  matrix->istep = (max - min)/(double)nbins;
  for(i=0;i < nbins;i++){
    *(matrix->matrixptr + i) = 0.0;
    *(matrix->imid + i) = min + matrix->istep/2.0 + (double)i*matrix->istep;
  }
  matrix->nibins = nbins;
  matrix->njbins = 0;
  matrix->imin = min; matrix->imax = max;

  return;
}

/****************************************************************************/
/* mbook2: set up 2-d histogram                                             */
/****************************************************************************/
void mbook2(MHist *matrix, int nibins, double imin, double imax, 
	    int njbins, double jmin, double jmax)
{
  int i, j;

  matrix->matrixptr = (double *)malloc(nibins*njbins*sizeof(double));
  matrix->imid = (double *)malloc(nibins*sizeof(double));
  matrix->jmid = (double *)malloc(njbins*sizeof(double));
  matrix->istep = (imax - imin)/(double)nibins;  
  matrix->jstep = (jmax - jmin)/(double)njbins;
  for(i=0;i < nibins;i++){
    *(matrix->imid + i) = imin + matrix->istep/2.0 + (double)i*matrix->istep;
    for(j=0;j < njbins;j++){
      if(i==0)
	*(matrix->jmid + j) = jmin + matrix->jstep/2.0 + 
	  (double)j*matrix->jstep;
      *(matrix->matrixptr + i + j*nibins) = 0.0;
    }
  }
  matrix->nibins = nibins; matrix->njbins = njbins;
  matrix->imin = imin;  matrix->jmin = jmin; 
  matrix->imax = imax;  matrix->jmax = jmax; 
  
  return;
}
