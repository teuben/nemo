/* -------------------------------------------------------------- *\
|* toolsnemo.c :	JC LAMBERT	le 24-Mar-95	V1.0
|*
|* Wrapper of basic NEMO procedure, it makes me life easiest :)
|*
|*
\* -------------------------------------------------------------- */

/* -------------------------------------------------------------- *\
|* Include files
\* -------------------------------------------------------------- */
#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>		
#include <filestruct.h>

#include <snapshot/snapshot.h>	
#include <snapshot/body.h>


/* -------------------------------------------------------------- *\ 
|* get_data_gen :
|* 
\* -------------------------------------------------------------- */
int get_data_gen(stream instr, char * TypeTag,char * DataType,
		 int size_alloc, int nbody, int dim1,int dim2,void ** genptr)
{
  if (*genptr == NULL)
    *genptr = ( void **) allocate(size_alloc);
  if (*genptr == NULL)
    error("Pas assez de memoire dans \"get_data_gen\"\n");
  get_data_coerced(instr, TypeTag, DataType,*genptr, nbody,
                   dim1,dim2,0);
}
/* -------------------------------------------------------------- *\ 
|* get_data_time :
|* Get Time from a Nemo snapshot
\* -------------------------------------------------------------- */
int get_data_time(stream instr, char * DataType, int size_type,
		  void ** timeptr)
{
  int etat;

  if (get_tag_ok(instr, TimeTag))
    { 
      if (*timeptr == NULL)
	*timeptr = (void **) allocate(size_type);
      get_data_coerced(instr, TimeTag, DataType, *timeptr,0);
      etat = 1; 
    }
  else
    etat = 0;
  return etat;
}
/* -------------------------------------------------------------- *\ 
|* get_data_nbody :
|* Get Nbody from aNEMO snapshot
\* -------------------------------------------------------------- */
int get_data_nbody(stream instr, char * DataType, int size_type, void ** nbodyptr)
{
  int etat;

  if (get_tag_ok(instr, NobjTag))
    { 
      if (*nbodyptr == NULL)
	*nbodyptr = (void **) allocate(size_type);
      get_data_coerced(instr, NobjTag, DataType, *nbodyptr,0);
      etat = 1; 
    }
  else
    etat = 0;
  return etat;
}
/* -------------------------------------------------------------- *\ 
|* get_data_mass :
|* Get Masses from a NEMO snapshot
\* -------------------------------------------------------------- */
int get_data_mass(stream instr, char * DataType, int nbody, int size_type,
		  void ** massptr)
{
  int etat;

  if (get_tag_ok(instr, MassTag))
    { 
      if (*massptr == NULL)
	*massptr = (void **) allocate(size_type * nbody);
      get_data_coerced(instr, MassTag, DataType, *massptr,
		       nbody, 0);
      etat = 1; 
    }
  else
    etat = 0;
  return etat;
}

/* -------------------------------------------------------------- *\ 
|* get_data_phase :
|* Get Phase Space coordinates from a NEMO snapshot
\* -------------------------------------------------------------- */
int get_data_phase(stream instr, char * DataType, int nbody, 
		   int size_type, void ** phaseptr, int ndim)
{
  int etat;

  if (get_tag_ok(instr, PhaseSpaceTag))
    { 
      if (*phaseptr == NULL)
	*phaseptr = (void **) allocate(size_type*nbody*2*ndim);

      get_data_coerced(instr, PhaseSpaceTag, DataType, *phaseptr,
		       nbody,2,ndim,0);

      etat = 1; 
    }
  else
    etat = 0;
  return etat;
}


/* -------------------------------------------------------------- *\ 
|* get_data_pot :
|* Get Potential from a NEMO snapshot
\* -------------------------------------------------------------- */
int get_data_pot(stream instr, char * DataType, int nbody, int size_type,
		 void ** potptr)
{
  int etat;

  if (get_tag_ok(instr, PotentialTag))
    { 
      if (*potptr == NULL)
	*potptr = (void **) allocate(size_type * nbody);
      get_data_coerced(instr, PotentialTag, DataType, *potptr,
		       nbody, 0);
      etat = 1; 
    }
  else
    etat = 0;
  return etat;
}
/* -------------------------------------------------------------- *\ 
|* get_data_acc :
|* Get Acceleration from a NEMO snapshot
\* -------------------------------------------------------------- */
int get_data_acc(stream instr, char * DataType, int nbody, int size_type,
		 void ** accptr, int ndim)
{
  int etat;

  if (get_tag_ok(instr, AccelerationTag))
    { 
      if (*accptr == NULL)
	*accptr = (void **) allocate(size_type * nbody*ndim);
      get_data_coerced(instr, AccelerationTag, DataType, *accptr,
		       nbody,ndim, 0);
      etat = 1; 
    }
  else
    etat = 0;
  return etat;
}
/* -------------------------------------------------------------- *\ 
|* chk_select : 
|*
|* return a double dimension array of boolean Q[nb_sel][nbody] 
|* 
|* 'nbody' is self explained
|* 'nb_sel' is the number of particle selection string (select_pts)
|* 'slect_pts' is an array of particle selection string like
|*              0:2000 || 0:2000:10 || 0:12000,23000:300000 etc...
|* 'nret' is an array of int of size nb_sel, which contain the number
|*        of selected particles for the current selection string
|*
|* 'Q' is a double dimension array of boolean which contain TRUE for
|*     the bodies index which  match to the selection string, otherwise
|      FALSE.
|*           
\* -------------------------------------------------------------- */
bool ** chk_select(int * nret,int nb_sel,int nbody,string select_pts[])
{
  int  ** select_i;
  bool ** Qsel;
  int i,j,k;

  /* Memory allocation */
  Qsel = (bool **) allocate(nb_sel*sizeof(int));
  select_i = (int **) allocate(nb_sel*sizeof(int));

  for (i=0; i<nb_sel; i++)
    {

      Qsel[i]  = (bool *) allocate(nbody*sizeof(bool));

      select_i[i] = (int * ) allocate(nbody*sizeof(int));
    }

  /* loop on all selection string */
  for (j=0; j<nb_sel; j++)
    {    

      /* Init Q array with FALSE */
      for (i=0 ; i < nbody ; i++)
	Qsel[j][i] = FALSE;

      if (!streq("all",select_pts[j]))
	{
	    
	  for (i=0; i<nbody; i++)
	    { 
	      Qsel[j][i] = FALSE;
	      select_i[j][i] =-1;
	    }
		
	  nret[j] = nemoinpi(select_pts[j],
			     select_i[j],nbody);

	  for (i=0; i < nret[j]; i++)
	    {
	      Qsel[j][select_i[j][i]]=TRUE;
	    }
	
	}
      else
	{ 
	  for (i=0; i<nbody; i++)
	    Qsel[j][i] = TRUE;   
	  nret[j] = nbody;
	}

    } 

  /* Free memory */
  for (k=0; k<nb_sel; k++)
    free((int **) select_i[k]);
  free((int *) select_i);


  return  Qsel;

}
/* -------------------------------------------------------------- *\ 
|* Fin de toolsnemo.c 
\* -------------------------------------------------------------- */ 
