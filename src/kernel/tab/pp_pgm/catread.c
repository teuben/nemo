#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "catread.h"
#include "catutil.h"

#define NSPC 512
static int cat_nmol = 0;
static struct {
  int moltag,nlen,ver;
  float qln[7];
  char cname[16];
} catcom[NSPC+1], *caterr, *catptr;

int catfrq(molec, cfreq, buff)
int molec;
char *cfreq, *buff;
{
  /*****************************************************************

    CATFRQ CALLS
    CATDIR CALLS CATRD
    CATRD CALLS CATFIL

    FUNCTION TO FIND AND RETURN FIRST LINE FOR SPECIES "MOLEC"
    THAT IS BEYOND FREQ

    CATFRQ <0 IF ERROR FOUND
    CATFRQ =0 IF FREQ GT FREQUENCY OF LAST LINE
    CATFRQ >0 POINTS TO LINE 

    BUFF CONTAINS THE LINE INFORMATION

    CFREQ IS CHARACTER STRING IN FORMAT (F13.4) OR EQUIVALENT

  ******************************************************************/
  int line,nline,l,r;
  nline = catlen(molec);   
  /*     lines  < L have frequency  < FREQ */
  /*     lines >= R have frequency >= FREQ */
  l = line = 1; r = nline + 1; buff[0] = 0;
  while(l < r) {
    line = (l + r) >>1;
    if (catrd(molec, line, buff)) return (-1);
    if (strncmp(cfreq, buff, 13) >= 0){
      l = line + 1;
    }else{
      r = line;
    }
  }
  if(l > nline) {
    buff[0] = 0; return (0);
  }
  if(r == line) return (line);
  line = r;
  if(catrd(molec,line,buff)) return (-1);
  return (line);
}

int catrd(molec, line, buff)
int molec, line;
char* buff;
/***********************************************************************

     CATRD CALLS
       CATFIL


     SUBROUTINE TO READ CATALOGUE FILE FOR MOLECULE "MOLEC",LINE# "LINE"

       80 CHARACTER BUFFER RETURNED IN BUFF

     ERROR CODE RETURNED (IOSTAT ERROR #)

*************************************************************************/
{
  static int buflen=80;
  static int omolec=-1;
  static FILE *handle;
  long offset;
  if (molec != omolec) { /* new molecule */
    if(omolec > 0) fclose(handle);
    omolec = molec; 
    handle = fopen(catfil(molec),"r");
    if(handle == NULL) {
      buff[0] = 0; omolec = -1; 
      return (-1);
    }
  }
  if(line <= 0) return(-2);
  offset = (line-1) * (long)(buflen + 1);
  fseek(handle, offset, SEEK_SET);
  if (fread(buff, 1, buflen, handle) == buflen) {
    buff[buflen]=0;
    return (0);
  }else{
    buff[0]=0;
    return (1);
  }
}

char *catfil(num)
int num;
/**********************************************************
     FUNCTION WHICH RETURNS FILE NAME CORRESPONDING TO NUM
***********************************************************/
{
  static char catdir[]={"/catalog/catdir.cat"};
  static char catent[]={"/catalog/c000000.cat"};
  char *cfield; 
  int k,iq;
  if(num == 0) return catdir;
  cfield = &catent[16];
  for (k = 0; k < 6; ++k){
    iq=num; num=num/10; --cfield;
    *cfield = (char)('0'+ (iq - num * 10));
  }
  return catent;
}
  
int catlen(molec)
int molec;
/**************************************************
C
C   SUBROUTINE TO RETURN CATALOGUE ENTRY LENGTH
C
C   MOLEC IS SPECIES TAG
C
****************************************************/
{
  static int fmt[9] = {6,7,7,7,7,7,7,7,2}; 
  FILE *cdir;
  char *pbuf;
  double dval[9];
  char buf[82];
  float *qptr;
  int k;
  if(cat_nmol == 0) { /* initialize */
    cdir = fopen(catfil(0),"r");
    if(cdir == NULL) return (0);
    catptr = catcom;
    while(cat_nmol < NSPC) {
      if(fgets(buf,82,cdir) == NULL) break;
      pcard(buf, dval, 1, fmt);
      catptr->moltag = (int)dval[0];  
      pbuf = catptr->cname;
      memcpy(pbuf, buf + 6, 14); pbuf[14] = 0;
      pcard(buf + 20, dval, 9, fmt);
      catptr->nlen = (int) dval[0];
      if (catptr->moltag == 0 || catptr->nlen == 0) continue;
      qptr = catptr->qln;
      for (k = 0; k < 7; k++) qptr[k] = (float) dval[k+1];
      catptr->ver = (int) dval[8]; 
      ++catptr; ++cat_nmol;      
    }
    fclose (cdir);
    caterr = catptr;
    strcpy(caterr->cname, " error");
    qptr = caterr->qln;
    for (k = 0; k < 7; k++) qptr[k] = 0.;
    caterr->moltag = 0; caterr->nlen = 0; caterr->ver = 0;
  }
  if(molec > 0) {
    for(k = 0; k < cat_nmol; ++k){
      if (catptr == caterr) catptr = catcom;
      if (molec == catptr->moltag) return catptr->nlen;
      ++catptr;
    }
  }
  catptr = caterr;
  return 0;
}

char *catdir(molec,nline,qqln,iver)
int molec, *nline, *iver;
double *qqln;
/*********************************************************
C   SUBROUTINE TO RETURN CATALOGUE DIRECTORY INFORMATION
C
C   MOLEC IS SPECIES TAG
C   CATDIR IS ASCII SPECIES NAME
C   NLINE IS THE NUMBER OF LINES FOR MOLECULE
C   QLN IS THE LOG10 OF THE PARTITION FUNCTION FOR
C        300,225,150,75,37.5,18.75,9.375 K
C   IVER IS THE VERSION NUMBER
*************************************************************/
{
  int k;
  float *qptr;
  *nline = catlen(molec);
  qptr=catptr->qln;
  for(k = 0; k < 7; k++) qqln[k] = qptr[k];
  *iver = catptr->ver;
  return catptr->cname;
}

int nxtdir(molec)
int *molec;
/**********************************************************************
C
C     FUNCTION NXTDIR RETURNS THE NUMBER OF REMAINING DIRECTORY ENTRIES
C        AFTER INPUT SPECIES MOLEC
C     ON RETURN MOLEC IS CHANGED TO THE NEXT SPECIES TAG
C
C     IF INPUT SPECIES MOLEC = 0 THEN POSITION IS BEFORE FIRST ENTRY
C
C     IF INPUT SPECIES MOLEC = LAST ENTRY THEN MOLEC=0 ON RETURN
***********************************************************************/
{
  catlen(*molec); 
  if (*molec == 0)
    catptr = catcom;
  else if (catptr != caterr) 
    ++catptr; 
  *molec = catptr->moltag;
  return catptr->nlen;
}

int getcat(buf, pdata)
char *buf;
struct catdata *pdata;
{
  static double dval[8];
  static int fmt[8]={13,8,8,2,10,3,7,4};
  if (pcard(buf, dval, 8, fmt) < 8) return -1;
  pdata->freq = dval[0];
  pdata->derr = dval[1];
  pdata->str  = dval[2];
  pdata->itd  = (int)dval[3];
  pdata->elow = dval[4];
  pdata->igup = (int)dval[5];
  pdata->tag  = (int)dval[6];
  pdata->ifmt = (int)dval[7];
  return (readqn(buf + 55, pdata->iqn, 12));
}
