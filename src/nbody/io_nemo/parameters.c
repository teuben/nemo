/* -------------------------------------------------------------- *\
|* parameters.c :	JC LAMBERT	 09-Feb-98	V1.1
|*
|* V1.0  : created
|*
|* V1.1  : minor modifications
|*
|* Materiel     : Sparc STATION
|* OS           : Solaris 2.x
|* Langage      : C
|* Version Nemo : 2.0 Juillet 1994 
|*
\* -------------------------------------------------------------- */

/* -------------------------------------------------------------- *\
|* Standard include files
\* -------------------------------------------------------------- */
#include <stdinc.h>
#include <string.h>

/* -------------------------------------------------------------- *\
|* Nemo include files
\* -------------------------------------------------------------- */
#include "flags_data.h"


/* -------------------------------------------------------------- *\ 
|* get_case :
|* Return the code of the fields if it exist otherwise 0.
\* -------------------------------------------------------------- */
int get_case(char * champs)
{
  struct ch_num
    { char * ch;
      int    num;
    } tab_ch_num[] =
        { "n",  '\001',  "nbody",   '\001',   /* nbody */
          "t",  '\002',  "time",    '\002',   /* time  */
          "m",  '\003',  "mass",    '\003',   /* mass  */
          "x",  '\004',  "pos",     '\004',   /* pos   */
          "v",  '\005',  "vel",     '\005',   /* vel   */
          "p",  '\006',  "pot",     '\006',   /* pot   */
          "close",  '\007',                   /* close */
          "s",  '\010',  "save",     '\010',  /* save  */
          "r",  '\011',  "read",     '\011',  /* read  */
          "a",  '\012',  "acc",      '\012',  /* acc   */
          "single", '\013', "float", '\013',  /* float type */
          "double", '\014', "double", '\014', /* double type */
	  "info",   '\015', "diag",   '\015', /* print information */
	  "h",      '\016', "history", '\016',/* history */
	  "st", '\017',                       /* selected time */
	  "sp", '\020',                     /* selected particles */
          NULL
        };   
  struct ch_num * p;
  int find = '\000';

  p = tab_ch_num;
  while(p->ch)
    { 
#ifdef DEBUG1
       fprintf(stderr,"[get_case] >%s< --- >%s< %d\n",champs,p->ch,
                       strcmp(p->ch,champs));
#endif
      if (!strcmp(p->ch,champs))
        { find = p->num;
#ifdef DEBUG1
       fprintf(stderr,"[get_case] cmp ok find = %d\n",find);
#endif
          break;
        }
      else
        p++;
    }

  return find;
}

/* -------------------------------------------------------------- *\ 
|* get_champs :
|* Return the field between two comma.
\* -------------------------------------------------------------- */
char * get_champs(char **p)
{ 
  char * chaine="", * x, * y;

  int len,start,end,i;

  /* find out the next ',' from the current position */
  x = strchr(*p,',');
  y = *p;
  if (x)
    len = x-*p;
  else
    len = strlen(*p);

  /* find out the boundaries of the string (without blank) */

  /* 'start' match to the beginning of the string */
  for (i=0; i < len; i++)
    {
      if (y[i] != ' ')
	{
	  start = i;
	  break;
	}
    }
  /* 'end' match to the end of the string */
  for (i=0; i < len; i++)
    {
      if (y[len-1-i] != ' ')
	{
	  end = len-1-i  ;
	  break;
	}
    }

  chaine = (char * ) malloc(sizeof(char)*(end-start+2));
  if (!chaine)
    { 
      fprintf(stderr,"[get_champs] memory allocation error.\n");
      exit(1);
    } 
  
  /* copy the selected string into chaine */
  strncpy(chaine,y+start,end-start+1);
  chaine[end-start+1]='\0';

  /* set the pointer to the next parameters */
   if (x)
    *p=x+1 ;
   else
    *p=*p+len;

  return(chaine);
}
/* -------------------------------------------------------------- *\ 
|* chk_parameters :
|* Print out selected parameters.
\* -------------------------------------------------------------- */
int chk_parameters(bool io_op, int rtype)
{ 

  char * tab_info_real[2] = { "Float", "Double"};

  if (io_op)
     fprintf(stderr,"Reading .... \n[");
  else
     fprintf(stderr,"Saving .... \n[");

  if (N_io)
     fprintf(stderr," n");
  if (T_io)
     fprintf(stderr," t");
  if (M_io)
     fprintf(stderr," m");
  if (X_io)
     fprintf(stderr," x");
  if (V_io)
     fprintf(stderr," v");
  if (P_io)
     fprintf(stderr," p");
  if (A_io)
     fprintf(stderr," a");
  fprintf(stderr," <%s> ]\n",tab_info_real[rtype-1]);

}
/* -------------------------------------------------------------- *\
|* char2float():
|* convert a string to a float
\* -------------------------------------------------------------- */
float char2float(char * ch,int rtype)
{
  float r;

  switch(rtype)
    {
    case 1 : r = (float) *((float *) ch);
      break;
 
    case 2 : r = (float) *((double *) ch);
      break;
    }
  return r;
}
/* -------------------------------------------------------------- *\
|* char2double():
|* convert a string to a double
\* -------------------------------------------------------------- */
double char2double(char * ch,int rtype)
{
  double r;

  switch(rtype)
    {
    case 1 : r = (double) *((float *) ch);
      break;
 
    case 2 : r = (double) *((double *) ch);
      break;
    }
  return r;
  
}
/* -------------------------------------------------------------- *\
|* End of parameters.c
\* -------------------------------------------------------------- */
