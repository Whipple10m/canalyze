#include <stdio.h>
#include <stdlib.h>
 
#include "hrc.h"

struct header get_head(FILE *infile)
{
  int itmp;

  struct header head;
  
  head.runid = (char *)malloc(sizeof(char)*4);
  head.mode = (char *)malloc(sizeof(char)*3);
  head.sourcename = (char *)malloc(sizeof(char)*20);
  head.skyq = (char *)malloc(sizeof(char)*6);
  head.comms = (char *)malloc(sizeof(char)*404);

  head.flag = 0;
  fseek(infile,4L,SEEK_SET);
  if((fread(head.runid,sizeof(char),4,infile) == 0)){
    fprintf(stderr,
	    "   **error** reading header! Probably not a reduced file.\n");
    head.flag = 1;
    return(head);
  }
  if((fread(&head.nevents,sizeof(int),1,infile) == 0)){
    fprintf(stderr,
	    "   **error** reading header! Probably not a reduced file.\n");
    head.flag = 1;
    return(head);
  }
  if((fread(&head.duration,sizeof(double),1,infile) == 0)){
    fprintf(stderr,
	    "   **error** reading header! Probably not a reduced file.\n");
    head.flag = 1;
    return(head);
  }
  vaxtoieee_double(&head.duration);
  if((fread(&head.kdur,sizeof(int),1,infile) == 0)){
    fprintf(stderr,
	    "   **error** reading header! Probably not a reduced file.\n");
    head.flag = 1;
    return(head);
  }
  if((fread(head.mode,sizeof(char),3,infile) == 0)){
    fprintf(stderr,
	    "   **error** reading header! Probably not a reduced file.\n");
    head.flag = 1;
    return(head);
  }
  if((fread(head.sourcename,sizeof(char),20,infile) == 0)){
    fprintf(stderr,
	    "   **error** reading header! Probably not a reduced file.\n");
    head.flag = 1;
    return(head);
  }
  if((fread(&head.kdate,sizeof(int),1,infile) == 0)){
    fprintf(stderr,
	    "   **error** reading header! Probably not a reduced file.\n");
    head.flag = 1;
    return(head);
  }
  if((fread(&head.mjd,sizeof(double),1,infile) == 0)){
    fprintf(stderr,
	    "   **error** reading header! Probably not a reduced file.\n");
    head.flag = 1;
    return(head);
  }
  vaxtoieee_double(&head.mjd);
  if((fread(&head.frjd,sizeof(double),1,infile) == 0)){
    fprintf(stderr,
	    "   **error** reading header! Probably not a reduced file.\n");
    head.flag = 1;
    return(head);
  }
  vaxtoieee_double(&head.frjd);
  if((fread(&head.ra,sizeof(float),1,infile) == 0)){
    fprintf(stderr,
	    "   **error** reading header! Probably not a reduced file.\n");
    head.flag = 1;
    return(head);
  }
  vaxtoieee_float(&head.ra);
  if((fread(&head.dec,sizeof(float),1,infile) == 0)){
    fprintf(stderr,
	    "   **error** reading header! Probably not a reduced file.\n");
    head.flag = 1;
    return(head);
  }
  vaxtoieee_float(&head.dec);
  if((fread(&head.kut,sizeof(int),1,infile) == 0)){
    fprintf(stderr,
	    "   **error** reading header! Probably not a reduced file.\n");
    head.flag = 1;
    return(head);
  }
  if((fread(&head.kst,sizeof(int),1,infile) == 0)){
    fprintf(stderr,
	    "   **error** reading header! Probably not a reduced file.\n");
    head.flag = 1;
    return(head);
  }
  if((fread(&head.azimuth,sizeof(float),1,infile) == 0)){
    fprintf(stderr,
	    "   **error** reading header! Probably not a reduced file.\n");
    head.flag = 1;
    return(head);
  }
  vaxtoieee_float(&head.azimuth);
  if((fread(&head.elevation,sizeof(float),1,infile) == 0)){
    fprintf(stderr,
	    "   **error** reading header! Probably not a reduced file.\n");
    head.flag = 1;
    return(head);
  }
  vaxtoieee_float(&head.elevation);
  if((fread(head.skyq,sizeof(char),6,infile) == 0)){
    fprintf(stderr,
	    "   **error** reading header! Probably not a reduced file.\n");
    head.flag = 1;
    return(head);
  }
  if((fread(head.comms,sizeof(char),404,infile) == 0)){
    fprintf(stderr,
	    "   **error** reading header! Probably not a reduced file.\n");
    head.flag = 1;
    return(head);
  }
  if((fread(&head.gpsbeg,sizeof(double),1,infile) == 0)){
    fprintf(stderr,
	    "   **error** reading header! Probably not a reduced file.\n");
    head.flag = 1;
    return(head);
  }

  return(head);
}

struct frame get_frame(FILE *infile){

  struct frame fr;

  fr.channel = (short *)malloc(sizeof(short)*124);

  fr.flag = 0;

  if(fread(&fr.code,sizeof(int),1,infile)==0){
    fr.flag = 1;
    return(fr);
  }
  if(fread(&fr.time,sizeof(double),1,infile)==0){
    fr.flag = 1;
    return(fr);
  }
  if(fread(fr.channel,sizeof(short),124,infile)==0)
    fr.flag = 1;
  vaxtoieee_double(&fr.time);
  return(fr);
}
