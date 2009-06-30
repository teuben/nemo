/* =================================================================
|  Copyright Jean-Charles LAMBERT - 2005                            
|  e-mail:   Jean-Charles.Lambert@oamp.fr                           
|  address:  Dynamique des galaxies                                 
|            Laboratoire d'Astrophysique de Marseille               
|            2, place Le Verrier                                    
|            13248 Marseille Cedex 4, France                        
|            CNRS U.M.R 6110                                        
| ==================================================================
|* Set of usefull functions                                         
+----------------------------------------------------------------- */
#include <stdinc.h>
#include <stdlib.h>
#include <getparam.h>

#include "io_nemo_tools.h"

/* ----------------------------------------------------------------
|  allocate_pointer :                                              
|  check if the pointer is already allocated. if not it is         
|  allocated.                                                      
+---------------------------------------------------------------- */ 
char * allocate_pointer(char * p, int lg)
{ char * ptr;

 if (p == NULL) {
   ptr = (char *) malloc(sizeof(char)*lg);
   if (!ptr) {
     fprintf(stderr,
	     "[allocate_pointer], allocation memory error, aborted\n");
     exit(1);
   }
   return ptr;
 }
 else /* assume p is already  allocated */
   return p;
} 
/* ----------------------------------------------------------------
|  char2float():                                                   
|  convert a string to a float                                     
+---------------------------------------------------------------- */
float char2float(char * ch,int rtype)
{
  float r;

  switch(rtype) {
  case 1 : r = (float) *((float *) ch);
    break;
 
  case 2 : r = (float) *((double *) ch);
    break;
  }
  return r;
}
/* ----------------------------------------------------------------
|  char2double():                                                  
|  convert a string to a double                                    
+---------------------------------------------------------------- */
double char2double(char * ch,int rtype)
{
  double r;

  switch(rtype) {
  case 1 : r = (double) *((float *) ch);
    break;
    
  case 2 : r = (double) *((double *) ch);
    break;
  }
  return r;
  
}

/* ----------------------------------------------------------------
|  chk_select :                                                    
|                                                                  
|  return a double dimension array of boolean Q[nb_sel][nbody]     
|                                                                  
|  'nbody' is self explained                                       
|  'nb_sel' is the number of particle selection string (select_pts)
|  'slect_pts' is an array of particle selection string like       
|               0:2000 || 0:2000:10 || 0:12000,23000:300000 etc... 
|  'nret' is an array of int of size nb_sel, which contain the     
|   number of selected particles for the current selection string  
|                                                                  
|  'Q' is a double dimension array of boolean which contain TRUE   
|      for the bodies index which  match to the selection string,  
|      otherwise FALSE.                                            
|                                                                  
+---------------------------------------------------------------- */
bool ** chk_select(int * nret,int nb_sel,int nbody,string select_pts[])
{
  int  ** select_i;
  bool ** Qsel;
  int i,j,k;
  
  /* Memory allocation */
  Qsel     = (bool **) allocate(nb_sel*sizeof(int));
  select_i = (int **) allocate(nb_sel*sizeof(int));

  for (i=0; i<nb_sel; i++) {
    Qsel[i]  = (bool *) allocate(nbody*sizeof(bool));
    select_i[i] = (int * ) allocate(nbody*sizeof(int));
  }

  /* loop on all selection string */
  for (j=0; j<nb_sel; j++) {
    for (i=0 ; i < nbody ; i++) {
      Qsel[j][i] = FALSE;                /* init Set Q=False */
    }

    if (!streq("all",select_pts[j])) {   /* if not all particles selected */
      for (i=0; i<nbody; i++) {
	Qsel[j][i] = FALSE;              /* Set Q=False for not selected particles */
	select_i[j][i] =-1;
      }
		
      nret[j] = nemoinpi(select_pts[j],
			 select_i[j],nbody);

      for (i=0; i < nret[j]; i++) {
	Qsel[j][select_i[j][i]]=TRUE;    /* Set Q=True for selected particles */
      }
    }
    else {
      for (i=0; i<nbody; i++) {          /* All the particles selected */
	Qsel[j][i] = TRUE;               /* Set Q=True for all particles */
      }  
      nret[j] = nbody;
    }
  } 

  /* Free memory */
  for (k=0; k<nb_sel; k++) {
    free((int **) select_i[k]);
  }
  free((int *) select_i);

  return  Qsel;

}
/* ----------------------------------------------------------------
|  f_ch_to_c :                                                     
|  Insert '\0' character to the 'lg' position of a FORTRAN string. 
+---------------------------------------------------------------- */ 
char * f_ch_to_c(char * chaine,int lg)
{
  char * tmp;
  
  /* chaine[lg] = '\0'; */
#if 1
  char *p;
  p = strchr(chaine,'\0');
  dprintf(1,"[f_ch_to_c] p=[%x] chaine=[%x] diff [%d] lg=<%d>\n",p,chaine,p-chaine,lg);
  if ((p - chaine) >= lg ) {
    dprintf(1,"[f_ch_to_c] gonna fix fortran supposed string...\n");
    tmp = chaine+lg-1;
    
    while (*tmp==' ') {
      *tmp='\0';
      tmp--;
    }
  }
#else
  tmp = chaine+lg-1;

  while (*tmp==' ') {
    *tmp='\0';
    tmp--;
  }  
#endif
  return(chaine);
} 
/* ----------------------------------------------------------------
|  get_selected :                                                  
|  Return the string limited by "#" character.                     
+---------------------------------------------------------------- */
char * get_selected(char *p)
{ 
  char * chaine=NULL, * x;
  int len;
  /* find out the '#' character */
  x = strchr(p,'#');
  if (!x) {
    fprintf(stderr,"[get_selected] error\n");
    fprintf(stderr,"You have forgotten to put a '#' at the end of a selected field (st or sp), aborted....\n");
    exit(1);
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
/* ----------------------------------------------------------------
|  set_eos :                                                  
|  set End Of String at the Character SEP
+---------------------------------------------------------------- */
char * set_eos(char *p, char sep)
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

/* ----------------------------------------------------------------
|  End of io_nemo_tools.c                                          
+---------------------------------------------------------------- */
