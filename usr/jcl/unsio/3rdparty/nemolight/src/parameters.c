/* =================================================================
|  Copyright Jean-Charles LAMBERT - 2008                            
|  e-mail:   Jean-Charles.Lambert@oamp.fr                           
|  address:  Dynamique des galaxies                                 
|            Laboratoire d'Astrophysique de Marseille               
|            2, place Le Verrier                                    
|            13248 Marseille Cedex 4, France                        
|            CNRS U.M.R 6110                                        
| ==================================================================
|* Manage io_nemo's parameters                                      
| ==================================================================
| 03-Mars-05: happy gcc 3.4                                         
+----------------------------------------------------------------- */

/* -----------------------------------------------------------------
|  Standard include files                                           
+----------------------------------------------------------------- */ 
#include <stdinc.h>

/* -----------------------------------------------------------------
|  Nemo include files                                               
+----------------------------------------------------------------- */
#include "parameters.h"

/* extern variables */
#include "flags_data.h"

/* -----------------------------------------------------------------
|  get_case :                                                       
|  Return the code of the fields if it exist otherwise 0.           
+----------------------------------------------------------------- */
int get_case(char * field)
{
  struct ch_num
    { char * ch;
      int    num;
    } tab_ch_num[] ={
    {"n"     ,   1}, {"nbody"    ,   1},   /* nbody              */
    {"t"     ,   2}, {"time"     ,   2},   /* time               */
    {"m"     ,   3}, {"mass"     ,   3},   /* mass               */
    {"x"     ,   4}, {"pos"      ,   4},   /* pos                */
    {"v"     ,   5}, {"vel"      ,   5},   /* vel                */
    {"p"     ,   6}, {"pot"      ,   6},   /* pot                */
    {"a"     ,   7}, {"acc"      ,   7},   /* acc                */
    {"k"     ,   8}, {"key"      ,   8},   /* keys               */
    {"xv"    ,   9}, {"phase"    ,   9},   /* phase              */
    {"e"     ,  10}, {"eps"      ,  10},   /* eps                */
    {"b"     ,  11}, {"bits"     ,  11},   /* bits               */
    {"X"     ,  12}, {"aux"      ,  12},   /* aux                */
    {"d"     ,  13}, {"dens"     ,  13},   /* dens               */
    {"f"     ,  50}, {"n3"       ,  50},   /* n3                 */
    {"c"     ,  51}, {"3n"       ,  51},   /* 3n                 */
    {"s"     ,  52}, {"save"     ,  52},   /* save               */
    {"r"     ,  53}, {"read"     ,  53},   /* read               */
    {"single",  54}, {"float"    ,  54},   /* float type         */
    {"real4" ,  54},                       /* float type         */
    {"double",  55}, {"double"   ,  55},   /* double type        */
    {"real8" ,  55},                       /* double type        */
    {"info"  ,  56}, {"diag"     ,  56},   /* print information  */
    {"st"    ,  57},                       /* selected time      */
    {"sp"    ,  58},                       /* Selected particles */
    {"h"     ,  59}, {"history"  ,  59},   /* history            */
    {"close" ,  60},                       /* close              */
    { NULL   ,  -1}
  };   
  struct ch_num * p;
  int find = '\000';

  p = tab_ch_num;
  while(p->ch) {
    if (!strcmp(p->ch,field)) {
      find = p->num;
      break;
    }
    else
      p++;
  }
  return find;
}

/* -----------------------------------------------------------------
|  get_field :                                                      
|  Return the field between two comma.                              
+----------------------------------------------------------------- */
char * get_field(char **p)
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
  for (i=0; i < len; i++) {
    if (y[i] != ' ') {
      start = i;
      break;
    }
  }
  /* 'end' match to the end of the string */
  for (i=0; i < len; i++) {
    if (y[len-1-i] != ' ') {
      end = len-1-i  ;
      break;
    }
  }

  chaine = (char * ) malloc(sizeof(char)*(end-start+2));
  if (!chaine) {
    fprintf(stderr,"[get_field] memory allocation error.\n");
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
/* -----------------------------------------------------------------
|  chk_parameters :                                                 
|  Print out selected parameters.                                   
+----------------------------------------------------------------- */
int chk_parameters(bool io_op,int size_array, int rtype)
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
  if (XV_io)
     fprintf(stderr," xv");
  if (P_io)
     fprintf(stderr," p");
  if (A_io)
     fprintf(stderr," a");
  if (AUX_io)
     fprintf(stderr," aux");
  if (EPS_io)
     fprintf(stderr," e");
  if (D_io)
     fprintf(stderr," d");
  if (K_io)
     fprintf(stderr," k");

  if ( size_array ) {  /* 'io_nemo_f' mode (FORTRAN) */
    if (F_dim)
      fprintf(stderr," Fortran(%d,3) <%s> ]\n",
	      size_array,tab_info_real[rtype-1]);
    else
      fprintf(stderr," Fortran(3,%d) <%s> ]\n",
	      size_array,tab_info_real[rtype-1]); 
  }
  else                  /* 'io_nemo' mode (C)          */
    fprintf(stderr," <%s> ]\n",tab_info_real[rtype-1]);

  return 1;
}

/* -----------------------------------------------------------------
|  End of parameters.c                                              
+----------------------------------------------------------------- */
