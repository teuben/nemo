/* ================================================================
|  Copyright Jean-Charles LAMBERT - 2008                           
|  e-mail:   Jean-Charles.Lambert@oamp.fr                          
|  address:  Dynamique des galaxies                                
|            Laboratoire d'Astrophysique de Marseille              
|            2, place Le Verrier                                   
|            13248 Marseille Cedex 4, France                       
|            CNRS U.M.R 6110                                       
| =================================================================
|* Wrapper of basic NEMO procedure, it makes me life easiest :)    
+---------------------------------------------------------------- */

/* ----------------------------------------------------------------
|  Include files                                                   
+---------------------------------------------------------------- */
#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>		
#include <filestruct.h>

#include <snapshot/snapshot.h>	
#include <snapshot/body.h>

extern int maxbodies[];
extern int CURRENT_IO;
/* ----------------------------------------------------------------
|  get_data_gen :                                                  
|                                                                  
+---------------------------------------------------------------- */
int get_data_gen(stream instr, char * TypeTag,char * DataType,
		 int size_alloc, int nbody, int dim1,int dim2,void ** genptr)
{
  if (*genptr && maxbodies[CURRENT_IO] < nbody) {
    free((void *) *genptr);
    *genptr = NULL;
  }
  if (*genptr == NULL)
    *genptr = ( void *) allocate(size_alloc);
  get_data_coerced(instr, TypeTag, DataType,*genptr, nbody,
                   dim1,dim2,0);
  return 0;
}
/* ----------------------------------------------------------------
|  get_data_time :                                                 
|  Get Time from a Nemo snapshot                                   
+---------------------------------------------------------------- */
int get_data_time(stream instr, char * DataType, int size_type,
		  void ** timeptr)
{
  int status;

  if (get_tag_ok(instr, TimeTag))
    { 
      if (*timeptr == NULL)
	*timeptr = (void *) allocate((size_t) size_type);
      get_data_coerced(instr, TimeTag, DataType, *timeptr,0);
      status = 1; 
    }
  else
    status = 0;
  return status;
}
/* ----------------------------------------------------------------
|  get_data_nbody :                                                
|  Get Nbody from aNEMO snapshot                                   
+---------------------------------------------------------------- */
int get_data_nbody(stream instr, char * DataType, int size_type, void ** nbodyptr)
{
  int status;

  if (get_tag_ok(instr, NobjTag))
    { 
      if (*nbodyptr == NULL)
	*nbodyptr = (void *) allocate((size_t) size_type);
      get_data_coerced(instr, NobjTag, DataType, *nbodyptr,0);
      status = 1; 
    }
  else
    status = 0;
  return status;
}
/* ----------------------------------------------------------------
|  get_data_mass :                                                 
|  Get Masses from a NEMO snapshot                                 
+---------------------------------------------------------------- */
int get_data_mass(stream instr, char * DataType, int nbody, int size_type,
		  void ** massptr)
{
  int status;

  if (get_tag_ok(instr, MassTag))
    { 
      if (*massptr && maxbodies[CURRENT_IO] < nbody) {
	dprintf(1,"NEW ALLOC => [%d] [%d]\n",maxbodies[CURRENT_IO],nbody);
	free((void *) *massptr);
	*massptr = NULL;
      }
      if (*massptr == NULL)
	*massptr = (void *) allocate((size_t) size_type * nbody);
      get_data_coerced(instr, MassTag, DataType, *massptr,
		       nbody, 0);
      status = 1; 
    }
  else
    status = 0;
  return status;
}
/* ----------------------------------------------------------------
|  get_data_pos :                                                  
|  Get Positions from a NEMO snapshot                              
+---------------------------------------------------------------- */
int get_data_pos(stream instr, char * DataType, int nbody, int size_type,
		  void ** posptr, int ndim)
{
  int status;

  if (get_tag_ok(instr, PosTag))
    { 
      if (*posptr && maxbodies[CURRENT_IO] < nbody) {
	dprintf(1,"pos NEW ALLOC => [%d] [%d]\n",maxbodies[CURRENT_IO],nbody);
	free((void *) *posptr);
	*posptr = NULL;
      }
      if (*posptr == NULL) {
	*posptr = (void *) allocate((size_t) size_type * nbody *ndim);
      }
      get_data_coerced(instr, PosTag, DataType, *posptr,
		       nbody,ndim,0);
      status = 1; 
    }
  else
    status = 0;
  return status;
}
/* ----------------------------------------------------------------
|  get_data_vel :                                                  
|  Get Velocities from a NEMO snapshot                             
+---------------------------------------------------------------- */
int get_data_vel(stream instr, char * DataType, int nbody, int size_type,
		  void ** velptr, int ndim)
{
  int status;

  if (get_tag_ok(instr, VelTag))
    { 
      if (*velptr && maxbodies[CURRENT_IO] < nbody) {
	free((void *) *velptr);
	*velptr = NULL;
      }
      if (*velptr == NULL)
	*velptr = (void *) allocate((size_t) size_type * nbody*ndim);
      get_data_coerced(instr, VelTag, DataType, *velptr,
		       nbody,ndim, 0);
      status = 1; 
    }
  else
    status = 0;
  return status;
}
/* ----------------------------------------------------------------
|  get_data_phase :                                                
|  Get Phase Space coordinates from a NEMO snapshot                
+---------------------------------------------------------------- */
int get_data_phase(stream instr, char * DataType, int nbody, 
		   int size_type, void ** phaseptr, int ndim)
{
  int status;

  if (get_tag_ok(instr, PhaseSpaceTag))
    { 
      if (*phaseptr && maxbodies[CURRENT_IO] < nbody) {
	free((void *) *phaseptr);
	*phaseptr = NULL;
      }
      if (*phaseptr == NULL)
	*phaseptr = (void *) allocate((size_t) size_type*nbody*2*ndim);

      get_data_coerced(instr, PhaseSpaceTag, DataType, *phaseptr,
		       nbody,2,ndim,0);

      status = 1; 
    }
  else
    status = 0;
  return status;
}


/* ----------------------------------------------------------------
|  get_data_pot :                                                  
|  Get Potential from a NEMO snapshot                              
+---------------------------------------------------------------- */
int get_data_pot(stream instr, char * DataType, int nbody, int size_type,
		 void ** potptr)
{
  int status;

  if (get_tag_ok(instr, PotentialTag))
    { 
      if (*potptr && maxbodies[CURRENT_IO] < nbody) {
	free((void *) *potptr);
	*potptr = NULL;
      }
      if (*potptr == NULL)
	*potptr = (void *) allocate((size_t) size_type * nbody);
      get_data_coerced(instr, PotentialTag, DataType, *potptr,
		       nbody, 0);
      status = 1; 
    }
  else
    status = 0;
  return status;
}
/* ----------------------------------------------------------------
|  get_data_acc :                                                  
|  Get Acceleration from a NEMO snapshot                           
+---------------------------------------------------------------- */
int get_data_acc(stream instr, char * DataType, int nbody, int size_type,
		 void ** accptr, int ndim)
{
  int status;

  if (get_tag_ok(instr, AccelerationTag))
    { 
      if (*accptr && maxbodies[CURRENT_IO] < nbody) {
	free((void *) *accptr);
	*accptr = NULL;
      }
      if (*accptr == NULL)
	*accptr = (void *) allocate((size_t) size_type * nbody*ndim);
      get_data_coerced(instr, AccelerationTag, DataType, *accptr,
		       nbody,ndim, 0);
      status = 1; 
    }
  else
    status = 0;
  return status;
}
/* ----------------------------------------------------------------
|  get_data_keys :                                                 
|  Get Keys from a NEMO snapshot                                   
+---------------------------------------------------------------- */
int get_data_keys(stream instr, char * DataType, int nbody, int size_type,
		  void ** keysptr)
{
  int status;

  if (get_tag_ok(instr, KeyTag))
    { 
      if (*keysptr && maxbodies[CURRENT_IO] < nbody) {
	free((void *) *keysptr);
	*keysptr = NULL;
      }
      if (*keysptr == NULL)
	*keysptr = (void *) allocate((size_t) size_type * nbody);
      get_data_coerced(instr, KeyTag, DataType, *keysptr,
		       nbody, 0);
      status = 1; 
    }
  else
    status = 0;
  return status;
}
/* ----------------------------------------------------------------
|  get_data_eps :                                                  
|  Get Eps from a NEMO snapshot                                    
+---------------------------------------------------------------- */
int get_data_eps(stream instr, char * DataType, int nbody, int size_type,
		  void ** epsptr)
{
  int status;

  if (get_tag_ok(instr, EpsTag))
    { 
      if (*epsptr && maxbodies[CURRENT_IO] < nbody) {
	free((void *) *epsptr);
	*epsptr = NULL;
      }
      if (*epsptr == NULL)
	*epsptr = (void *) allocate((size_t) size_type * nbody);
      get_data_coerced(instr, EpsTag, DataType, *epsptr,
		       nbody, 0);
      status = 1; 
    }
  else
    status = 0;
  return status;
}
/* ----------------------------------------------------------------
|  get_data_aux :                                                  
|  Get Aux from a NEMO snapshot                                    
+---------------------------------------------------------------- */
int get_data_aux(stream instr, char * DataType, int nbody, int size_type,
		  void ** auxptr)
{
  int status;

  if (get_tag_ok(instr, AuxTag))
    { 
      if (*auxptr && maxbodies[CURRENT_IO] < nbody) {
	free((void *) *auxptr);
	*auxptr = NULL;
      }
      if (*auxptr == NULL)
	*auxptr = (void *) allocate((size_t) size_type * nbody);
      get_data_coerced(instr, AuxTag, DataType, *auxptr,
		       nbody, 0);
      status = 1; 
    }
  else
    status = 0;
  return status;
}
/* ----------------------------------------------------------------
|  get_data_dens :                                                 
|  Get Dens from a NEMO snapshot                                   
+---------------------------------------------------------------- */
int get_data_dens(stream instr, char * DataType, int nbody, int size_type,
		  void ** densptr)
{
  int status;

  if (get_tag_ok(instr, DensityTag))
    { 
      if (*densptr && maxbodies[CURRENT_IO] < nbody) {
	free((void *) *densptr);
	*densptr = NULL;
      }
      if (*densptr == NULL)
	*densptr = (void *) allocate((size_t) size_type * nbody);
      get_data_coerced(instr, DensityTag, DataType, *densptr,
		       nbody, 0);
      status = 1; 
    }
  else
    status = 0;
  return status;
}
/* ----------------------------------------------------------------
|  End of [get_dat_nemo.c]                                         
+---------------------------------------------------------------- */ 
