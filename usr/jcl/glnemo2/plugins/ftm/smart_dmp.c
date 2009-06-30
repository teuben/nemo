/* -------------------------------------------------------------- *\
|* smart_dmp.c                                   JCL: 13-Feb-2001
|*
|* Read and convert DMP file to current binary architecture
|* (Convert little_endian <-> big_endian )
|*
\* -------------------------------------------------------------- */
#include <stdio.h>
#include <stdlib.h>

static int aswap=1;       /* auto swap */
static int detect_header=1;
static int fboundary; /* Fortran extra character */
static FILE * dmp_in;     /* DMP file handler */

#define PRINT_CHAR(x) fprintf(stderr,#x" [%s]\n",x)
#define PRINT_INT(x) fprintf(stderr,#x" [%d]\n",x)
#define PRINT_DOUBLE(x) fprintf(stderr,#x" [%f]\n",x)

/* Read 4 bytes character which surround Fortran unformated record */
#define READ_FMARK(x) fread(&fboundary,sizeof(int),1,x)
char *  set_eos2(char *,char);
/* -------------------------------------------------------------- *\
|* b_swap :
|* reverse byte order for the word 'x' 
|* example : x a 4 bytes word
|*    x     -> 1 2 3 4
|* b_swap x -> 4 3 2 1
\* -------------------------------------------------------------- */
void b_swap( void * x,
             int size /* size of the word */
                 )
{
  char * p, t;
  int i;
  p = x;
  for (i=0;i<size/2;i++)
   { 
     t=*(p+i);
     *(p+i)=*(p+size-i-1);
     *(p+size-i-1)=t;
   }
}
/* -------------------------------------------------------------- *\
|* smart_fread :
|* if auto_swap is activated this function will swap data from/to
|* littleendian/bigendian
\* -------------------------------------------------------------- */
size_t smart_fread(
          char *pp,
          size_t size,
          size_t num_items,
          FILE *stream)
{
  int i,is,*p;
  static int first=1;
  int fboundary;

  fread(pp,size,num_items,stream);

  if (detect_header && aswap)
    {
      /* try to check in reverse order */
      p = (int *) pp;
      is = *p;
      b_swap((int *) &is,sizeof(int));
      if (is > *p) /* hummm... seems NOT to be good in reverse order */
	{
	  aswap = 0; 
	}
      else /* SWAP is necessary */
	{
	  if (first)
	    fprintf(stderr,"\n\n-*- Auto swapping activated -*-\n\n");
	  
	  first = 1;
	}

    }

  if (aswap && (size > 1))
    {
      for (i=0; i<num_items; i++)
	{
	  b_swap(&pp[i*size],size);
	}
    }
  if (size == 1 ) /* assume that it's a FORTRAN string ... */
    {
      /* ...must put a '\0' at the end of the string */
    }
  
}
/* -------------------------------------------------------------- *\
|* read_dmp_header:
|* Read first DMP record and detect binary format
\* -------------------------------------------------------------- */
int read_dmp_header_(char   * dmp_name,
		     int    * dmpindx,
		     char   * dmp_date,
		     int    * ndim,
		     int    * eqnindx,
		     double * version,
		     double * gamma,
		     double * poly,
		     int    * lsfm,
		     double * tframe,
		     int    * n1,
		     int    * n2)
{

  /* Open DMP file */
  dmp_name = set_eos2(dmp_name,'\\');
  if ( ! (dmp_in = fopen(dmp_name,"rb")))
    {
      fprintf(stderr,"Unable to read DMP file [%s], aborted...\n",
	      dmp_name);
      exit(1);
    }

  READ_FMARK(dmp_in); /* < */
  detect_header = 1;
  smart_fread((char *) dmpindx,sizeof(int),1,dmp_in);
  detect_header = 0;
  smart_fread(dmp_date,sizeof(char),24,dmp_in);
  smart_fread((char *) ndim,sizeof(int),1,dmp_in);
  smart_fread((char *) eqnindx,sizeof(int),1,dmp_in);
  smart_fread((char *) version,sizeof(double),1,dmp_in);
  smart_fread((char *) gamma,sizeof(double),1,dmp_in);
  smart_fread((char *) poly,sizeof(double),1,dmp_in);
  smart_fread((char *) lsfm,sizeof(int),1,dmp_in);
  READ_FMARK(dmp_in); /* > */

  READ_FMARK(dmp_in); /* < */
  smart_fread((char *) tframe,sizeof(double),1,dmp_in);
  smart_fread((char *) n1,sizeof(int),1,dmp_in);
  smart_fread((char *) n2,sizeof(int),1,dmp_in);
  READ_FMARK(dmp_in); /* > */


}
/* -------------------------------------------------------------- *\
|* read_dmp_anybyte:
|* read array of 'n' dmpof size_byte 'data'
\* -------------------------------------------------------------- */
int read_dmp_anybyte_(void * data, int * n,int * size_byte, int  * skip)
{
  int status;
  if ((*skip)) {
    /*fprintf(stderr,"< ftell %d\n",ftell(dmp_in));*/
    status = fseek(dmp_in,(*size_byte)*(*n)+8,SEEK_CUR);
    /*fprintf(stderr,"> ftell %d\n",ftell(dmp_in));*/
  } else {
    READ_FMARK(dmp_in); /* < */
    smart_fread(data,*size_byte,*n,dmp_in);
    READ_FMARK(dmp_in); /* > */
  }
}
/* -------------------------------------------------------------- *\
|* read_dmp_close:
|* Close DMP file
\* -------------------------------------------------------------- */
int read_dmp_close_()
{
  fclose(dmp_in);
}
/* ----------------------------------------------------------------
|  set_eos :                                                  
|  set End Of String at the Character SEP
+---------------------------------------------------------------- */
char * set_eos2(char *p, char sep)
{ 
  char * chaine=NULL, * x;
  int len;
  /* find out the '\' character */
  x = strchr(p,sep);
  if (!x) {
    return p;
  }
  /* size of the real field */
  len = x-p+1;
	
  /* copy the real field into chaine */
  if (len) {
    chaine = (char * ) allocate(sizeof(char) * len + 1);
    strncpy(chaine,p,len-1);
    chaine[len-1] = '\0';
    
    /*fprintf(stderr,"[%s] selected field : [%s] \n",p,chaine);*/
  }
	
  return chaine;
}
 
/* -------------------------------------------------------------- *\
|* Main program 
\* -------------------------------------------------------------- */
#ifdef TEST
void main(int argc, char ** argv)
{ 
  char    dmp_name[50];
  int     dmpindx;
  char    dmp_date[50];
  int     ndim;
  int     eqnindx;
  double  version;
  double  gamma;
  double  poly;
  int     lsm;
  double  tframe;
  int     n1;
  int     n2;

  if (argc != 2)
    {
      fprintf(stderr,"Usage : %s dmp_file\n",argv[0]);
      exit(1);
    }

 read_dmp_header_(argv[1],&dmpindx,dmp_date,&ndim,
		  &eqnindx,&version,&gamma,&poly,
		  &lsm,&tframe,&n1,&n2);
}
#endif
/* -------------------------------------------------------------- *\
|* End of [smart_dmp.c]
\* -------------------------------------------------------------- */
