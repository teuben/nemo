#define CATDATA struct catdata
CATDATA {
  double freq, derr, str, elow;
  int itd, igup, tag, ifmt;
  short iqn[12];
};

int catfrq(int molec,char *cfreq, char *buff);
/*****************************************************************
C   FUNCTION TO FIND AND RETURN FIRST LINE FOR SPECIES "MOLEC"
C       THAT IS BEYOND FREQ
C
C   CATFRQ <0 IF ERROR FOUND
C   CATFRQ =0 IF FREQ GT FREQUENCY OF LAST LINE
C   CATFRQ >0 POINTS TO LINE
C
C   BUFF CONTAINS THE LINE INFORMATION
C
C   FREQ IS CHARACTER STRING IN FORMAT (F13.4) OR EQUIVALENT
C
******************************************************************/

int catrd(int molec,int line,char *buff);
/************************************************************************
C     SUBROUTINE TO READ CATALOG FILE FOR MOLECULE "MOLEC",LINE# "LINE"
C
C       80 CHARACTER BUFFER RETURNED IN BUFF
C
C     ERROR CODE RETURNED (IOSTAT ERROR #)
C
*************************************************************************/

char *catfil(int num);
/**********************************************************
     FUNCTION WHICH RETURNS FILE NAME CORRESPONDING TO NUM
***********************************************************/

int catlen(int molec);
/**************************************************
C
C   SUBROUTINE TO RETURN CATALOGUE ENTRY LENGTH
C
C   MOLEC IS SPECIES TAG
C
****************************************************/

char *catdir(int molec,int *nline,double *qqln,int *iver);
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

int nxtdir(int *molec);
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
int getcat(char *buf, struct catdata *pdata);
/* fill structure with catlog data */
