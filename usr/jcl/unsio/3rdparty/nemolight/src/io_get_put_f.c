/* =================================================================
|  Copyright Jean-Charles LAMBERT - 2008                            
|  e-mail:   Jean-Charles.Lambert@oamp.fr                           
|  address:  Dynamique des galaxies                                 
|            Laboratoire d'Astrophysique de Marseille               
|            2, place Le Verrier                                    
|            13248 Marseille Cedex 4, France                        
|            CNRS U.M.R 6110                                        
| ==================================================================
|* Get/Put selected Data from io_nemo_f                             
| ==================================================================
+----------------------------------------------------------------- */

/* ----------------------------------------------------------------
|  Nemo include files                                              
+---------------------------------------------------------------- */
#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>		
#include <filestruct.h>
#include <history.h>

#include <snapshot/snapshot.h>	
#include <snapshot/body.h>

/* ----------------------------------------------------------------
|  Local include files                                             
+---------------------------------------------------------------- */
#include "history_data_f.h"
#include "flags_data.h"
#include "io_nemo_data_f.h"

#include "io_get_put_f.h"
#include "check_file.h"
#include "io_nemo_tools.h"
#include "get_data_wrapper.h"
#include "parameters.h"


#define TIMEFUZZ 0.0000001
extern int maxbodies[];
int CURRENT_IO;
/* ----------------------------------------------------------------
|  put_data_select :                                               
|  Save the snapshot in NEMO format.                               
+---------------------------------------------------------------- */
int put_data_select_f(char * outfile,
		      int  * size_array,
		      int    rtype,
		      char * io_out[],
		      bool * save_one,
		      FILE * outstr[],
		      int    MAXIO)
{
  int coordsys =	CSCode(Cartesian, NDIM, 2),
    i,k,no_io;
  /* char * phasep;*/	/* Phase space coordinates */
  char * accptr;	/* Accelerations array     */
  char * posptr;	/* Positions array         */
  char * velptr;	/* Velocities array        */
	
  int jump = rtype*4,jump3=3*jump; /*jump6=6*jump;*/
  char *  OutType;
	
  if (rtype==1)
    OutType=FloatType;
  else
    OutType=DoubleType;
	
    
  /* Open the file for writing */
  if ((no_io = get_old_file(outfile,io_out,save_one,outstr,MAXIO)) < 0)
    no_io = get_new_file(outfile,io_out,save_one,outstr,"w",MAXIO);
	
  /* print out what is doing */
  if (I_io)
    chk_parameters(FALSE,*size_array, rtype);
	
#if 0	
  if (XV_io) {
    /* memory allocation for positions and velocities */
    phasep = (char *) malloc(sizeof(char)*2*NDIM*(*nbody_f)*jump);
    
    if (phasep == NULL) {
      fprintf(stderr,"Memory error ## [put_data_select]\n");
      fprintf(stderr,"Impossible to allocate PhaseSpace Coordinates\n");
      exit(1);
    }
		
    /* copy positions and velocities from FORTRAN to NEMO */
    if (!F_dim) { /* Fortran array : Tab(3,MAXBODY) */
      for(i=0; i<*nbody_f; i++) {
	if (X_io)
	  memcpy((char *) (phasep+i*jump6),
		 (char *) (pos_f+i*3*jump),
		 jump3);
				
	if (V_io)
	  memcpy((char *) (phasep+i*jump6+jump3),
		 (char *) (vel_f+i*3*jump),
		 jump3);
      }						  
    }
    else {		 /* Fortran array : Tab(MAXBODY,3) */
      for(i=0; i<*nbody_f; i++) {
	for(k=0; k<3; k++) {
	  if (X_io)
	    memcpy((char *) (phasep+i*jump6+k*jump),
		   (char *) (pos_f+k*(*size_array)*jump+i*jump),
		   jump);
	    
	  if (V_io)
	    memcpy((char *) (phasep+i*jump6+jump3+k*jump),
		   (char *) (vel_f+k*(*size_array)*jump+i*jump),
		   jump);
	}
      }
    }
  }
#endif
   
  if (X_io) {
    /* memory allocation for positions array */
    posptr = (char *) malloc(sizeof(char)*NDIM*(*nbody_f)*jump);
    if (posptr == NULL) {
      fprintf(stderr,"Memory error ## [put_data_select]\n");
      fprintf(stderr,"Impossible to allocate Posisitions array\n");
    }
		
    /* transfer acceleration */
    if (!F_dim) { /* Fortran array : Tab(3,MAXBODY) */
      for(i=0; i<*nbody_f; i++) {
	memcpy((char *)(posptr+i*3*jump),
	       (char *)(pos_f+i*3*jump),
	       jump3);
      } 						
    }
    else {		 /* Fortran array : Tab(MAXBODY,3) */
      for(i=0; i<*nbody_f; i++) {
	for(k=0; k<3; k++) {
	  memcpy((char *) (posptr+i*3*jump+k*jump),
		 (char *) (pos_f+k*(*size_array)*jump+i*jump),
		 jump);
	}
      }
    }
  }
	
  if (V_io) {
    /* memory allocation for velocities array */
    velptr = (char *) malloc(sizeof(char)*NDIM*(*nbody_f)*jump);
    if (velptr == NULL) {
      fprintf(stderr,"Memory error ## [put_data_select]\n");
      fprintf(stderr,"Impossible to allocate Velocities array\n");
    }
		
    /* transfer acceleration */
    if (!F_dim) { /* Fortran array : Tab(3,MAXBODY) */
      for(i=0; i<*nbody_f; i++) {
	memcpy((char *)(velptr+i*3*jump),
	       (char *)(vel_f+i*3*jump),
	       jump3);
      } 						
    }
    else {		 /* Fortran array : Tab(MAXBODY,3) */
      for(i=0; i<*nbody_f; i++) {
	for(k=0; k<3; k++) {
	  memcpy((char *) (velptr+i*3*jump+k*jump),
		 (char *) (vel_f+k*(*size_array)*jump+i*jump),
		 jump);
	}
      }
    }
  }

  /* memory allocation for acceleration */
  if (A_io) {
    accptr = (char *) malloc(sizeof(char)*NDIM*(*nbody_f)*jump);
    if (accptr == NULL) {
      fprintf(stderr,"Memory error ## [put_data_select]\n");
      fprintf(stderr,"Impossible to allocate Acceleration array\n");
    }
		
    /* transfer acceleration */
    if (!F_dim) { /* Fortran array : Tab(3,MAXBODY) */
      for(i=0; i<*nbody_f; i++) {
	memcpy((char *)(accptr+i*3*jump),
	       (char *)(acc_f+i*3*jump),
	       jump3);
      } 						
    }
    else {		 /* Fortran array : Tab(MAXBODY,3) */
      for(i=0; i<*nbody_f; i++) {
	for(k=0; k<3; k++) {
	  memcpy((char *) (accptr+i*3*jump+k*jump),
		 (char *) (acc_f+k*(*size_array)*jump+i*jump),
		 jump);
	}
      }
    }
  }
	
  /* save the snapshot */
  put_set(outstr[no_io], SnapShotTag);	
	
  put_set(outstr[no_io], ParametersTag);
  if (T_io)
    put_data(outstr[no_io], TimeTag, OutType, timu_f, 0);
  put_data(outstr[no_io], NobjTag, IntType, nbody_f, 0);
  put_tes(outstr[no_io], ParametersTag);
	   
  put_set(outstr[no_io], ParticlesTag); /*   start particle output  */
  put_data(outstr[no_io], CoordSystemTag, IntType, &coordsys, 0);

  if (M_io)
    put_data(outstr[no_io], MassTag, OutType, mass_f, *nbody_f, 0);

#if 0
  if (XV_io) {
    put_data(outstr[no_io], PhaseSpaceTag,OutType,phasep, *nbody_f, 
	     2, NDIM, 0); 
    free((char *) phasep); 
  }
#endif
  if (X_io) {
    put_data(outstr[no_io], PosTag , OutType, posptr, *nbody_f,
	     NDIM, 0);
    free((char *) posptr);
  }

  if (V_io) {
    put_data(outstr[no_io], VelTag , OutType, velptr, *nbody_f,
	     NDIM, 0);
    free((char *) velptr);
  }
  
  if (P_io)
    put_data(outstr[no_io], PotentialTag, OutType, pot_f, *nbody_f, 0);
  
  if (A_io) {
    put_data(outstr[no_io], AccelerationTag , OutType, accptr, *nbody_f,
	     NDIM, 0);
    free((char *) accptr);
  }
  if (AUX_io)
    put_data(outstr[no_io], AuxTag, OutType, aux_f, *nbody_f, 0);

  if (K_io)
    put_data(outstr[no_io], KeyTag, IntType, keys_f, *nbody_f, 0);

  if (D_io)
    put_data(outstr[no_io], DensityTag, OutType, dens_f, *nbody_f, 0);

  if (EPS_io)
    put_data(outstr[no_io], EpsTag, OutType, eps_f, *nbody_f, 0);

  put_tes(outstr[no_io], ParticlesTag);
	   
  put_tes(outstr[no_io], SnapShotTag);
	   
  save_one[no_io] = TRUE;
  return 1;
}
/* ----------------------------------------------------------------
|  get_data_select :                                               
|  Read the snapshot in NEMO format.                               
+---------------------------------------------------------------- */
int get_data_select_f(char * infile,
		      int  * size_array,
		      int    rtype,
		      char * io_in[],
		      bool * read_one,
		      FILE * instr[],
		      int    MAXIO) 
{
  char 
    * timeptr  = NULL,
    * phaseptr = NULL,
    * posptr   = NULL,
    * velptr   = NULL,
    * potptr   = NULL,
    * accptr   = NULL,
    * massptr  = NULL,
    * keysptr  = NULL,
    * auxptr   = NULL,
    * densptr  = NULL,
    * epsptr   = NULL,
    * SelString= NULL;
  int 
    * nbodyptr = NULL;	
	
  int * SelectedPart = NULL; /* added, table of indexes */
  int nBodySelected  = 0;
  int nParticlesRead = 0;
  char * select_time = NULL;

  double real_time;

  int i,k, no_io,i_jump,jump=rtype*4,jump3=3*jump,jump6=6*jump;
  char * OutType;
	
  string headline;

  i_jump = sizeof(int);  /* sizeof an Integer Jump */

  /* select the real format */
  if (rtype==1)
    OutType=FloatType;
  else
    OutType=DoubleType;
	
  /* Open the file for reading */
  if ((no_io = get_old_file(infile,io_in,read_one,instr,MAXIO)) < 0) {
    no_io = get_new_file(infile,io_in,read_one,instr,"r",MAXIO);
  }
  CURRENT_IO = no_io;

  /* print out what it's doing */
  if (I_io)
    chk_parameters(TRUE,*size_array,rtype);
	
  /* check the selected time field */
  if (ST_io) {
    select_time = (char *) get_selected(selt_f);
  }
	
  /* read the snapshot */
  for (;;) {
    get_history(instr[no_io]);	

    while (get_tag_ok(instr[no_io],HeadlineTag))
      headline = get_string(instr[no_io],HeadlineTag);
		
    if (!get_tag_ok(instr[no_io], SnapShotTag)) {
      if (!read_one[no_io]) { /* file has never been read */
	fprintf(stderr,"SnapshotTag error ## [get_data_select]\n");
	fprintf(stderr,"%s is not a NEMO SNAPSHOT\n",io_in[no_io]);
	exit(1);
      }

      /* end of snapshot reached */
      fprintf(stderr,"WARNING!! end of snapshot reached.\n");
      
      return 0;
    }
    else
      read_one[no_io] = TRUE; /* file has been read at least one time */
		
    get_set(instr[no_io], SnapShotTag);
		
    get_set(instr[no_io], ParametersTag);

    if (T_io) {  /* snapshot time */
      get_data_time(instr[no_io], OutType, rtype*4,
		    &timeptr);
      memcpy((char *) timu_f,(char *) timeptr, jump);
    }

    get_data_nbody(instr[no_io], IntType, sizeof(int),
		   &nbodyptr);
    if (SP_io) {
      /*get corrected selection string */
      SelString=(char *) get_selected(selp_f);

      if ( strcmp(SelString,"all")) {
	SelectedPart=(int*)allocate(*nbodyptr*sizeof(int));

	/*create index list */
	nBodySelected=nemoinpi(SelString,SelectedPart,*nbodyptr);
	if (nBodySelected < 0) {
	  fprintf(stderr,"Failed to select particles's range <%s>"
		  " *nemoinpi* function return code = [%d],"
		  " aborted.....\n",SelString,nBodySelected);
	  exit(1);
	}
#if 0
	fprintf(stderr,"select <%s>, n <%d>\n",
		SelString,nBodySelected);
#endif	
      }
      else { /* all particles are selected, so we *DESACTVATE* SP_io */
	SP_io = 0;
      }
    }
      
    if (N_io) /* nbody */
      *nbody_f = *nbodyptr;

		
    if (SP_io)
      nParticlesRead = nBodySelected;
    else
      nParticlesRead = *nbodyptr;

    /* test if FORTRAN array is enough big to get the NEMO data */
    if (*size_array < nParticlesRead) {
      fprintf(stderr,"Array error ## [get_data_select]\n");
      fprintf(stderr,"Array Fortran < Number of particles \n");
      fprintf(stderr,"Fortran max=(%d) <----> #particles read=(%d)\n",
	      *size_array,nParticlesRead);
      exit(1);				 
    }
		
    if (ST_io) { /* selected time */
      /* get the real value of the selected time */
      real_time= char2double(timu_f,rtype);
      
      if (!streq(select_time,"all") &&
	  !within(real_time,select_time,TIMEFUZZ)) {
	dprintf(1,"Info : skipping time step [%.4f]\n",real_time);
	get_tes(instr[no_io], ParametersTag);
	get_tes(instr[no_io], SnapShotTag);
	if (SP_io) {
	  free ((int *) SelectedPart);
	  free ((char *) SelString);
	}
	continue; /* go to the next time step */
      }
    }
		
    get_tes(instr[no_io], ParametersTag);
		
    if (get_tag_ok(instr[no_io], ParticlesTag)) {
      get_set(instr[no_io], ParticlesTag);
            
      /* get masses */
      if (M_io) {
	dprintf(1,"Getting Masses....\n");
	if (!get_data_mass(instr[no_io],OutType,*nbodyptr,
			   rtype*4,&massptr)) {
	  fprintf(stderr,"Snap error ### No Mass\n");
	  exit(1);
	} 	
	else {
	  if (SP_io)
	    for (i=0; i<nBodySelected; i++) {
	      memcpy((char*)(mass_f+jump*i),
		     (char*)(massptr+(jump*SelectedPart[i])),
		     jump);
	    }
	  else
	    memcpy((char *) (mass_f),
		   (char *) (massptr),
		   *nbodyptr * jump);
	  /*fprintf(stderr,"\n\nmass adress[%x]\n\n",massptr);*/
	  free((char *) massptr);
	}
      }  
      /* pos or/and vel selected */
      if (X_io || V_io  ) {
	dprintf(1,"Getting Phase Space....\n");
	if (get_data_phase(instr[no_io],OutType,*nbodyptr,
			   rtype*4,&phaseptr,NDIM)) {

	  if (SP_io) { /* particles selected */
	    if (!F_dim) { /* Fortran array : Tab(3,MAXBODY) */
	      for(i=0; i<nBodySelected; i++) {
		if (X_io) 
		  memcpy((char *) (pos_f+i*3*jump),
			 (char *) (phaseptr+SelectedPart[i]*jump6),
			 jump3);
		
		if (V_io)
		  memcpy((char *) (vel_f+i*3*jump),
			 (char *) (phaseptr+SelectedPart[i]*jump6+jump3),
			 jump3);
	      } 						
	    }
	    else {	   /* Fortran array : Tab(MAXBODY,3) */
	      for(i=0; i<nBodySelected; i++) {
		for(k=0; k<3; k++) {
		  if (X_io)
		    memcpy((char *) (pos_f+k*(*size_array)*jump+i*jump),
			   (char *) (phaseptr+SelectedPart[i]*jump6+k*jump),
			   jump);			   
										
		  if (V_io)
		    memcpy((char *) (vel_f+k*(*size_array)*jump+i*jump),
			   (char *) (phaseptr+SelectedPart[i]*jump6+jump3+k*jump),
			   jump); 
		}
	      }
	    }
	  } /* get_data_phase && SP_io */
	  else { /* get_data_phase && !SP_io */

	    /* copy positions and velocities from NEMO to FORTRAN */
	    if (!F_dim) { /* Fortran array : Tab(3,MAXBODY) */
	      for(i=0; i<*nbodyptr; i++) {
		if (X_io) 
		  memcpy((char *) (pos_f+i*3*jump),
			 (char *) (phaseptr+i*jump6),
			 jump3);

		
		if (V_io)
		  memcpy((char *) (vel_f+i*3*jump),
			 (char *) (phaseptr+i*jump6+jump3),
			 jump3);							
	      } 						
	    }
	    else {	   /* Fortran array : Tab(MAXBODY,3) */
	      for(i=0; i<*nbodyptr; i++) {
		for(k=0; k<3; k++) {
		  if (X_io)
		    memcpy((char *) (pos_f+k*(*size_array)*jump+i*jump),
			   (char *) (phaseptr+i*jump6+k*jump),
			   jump);			   
		  
		  if (V_io)
		    memcpy((char *) (vel_f+k*(*size_array)*jump+i*jump),
			   (char *) (phaseptr+i*jump6+jump3+k*jump),
			   jump); 
		}
	      }
	    }
	  } /* get_data_phase && SP_io */
	  free(phaseptr);
	} /* if get_data_phase */
	else { /* get_data_pos || get_data_vel */

	  if (X_io) {
	    dprintf(1,"Getting Positions....\n");
	    if (!get_data_pos(instr[no_io],OutType,*nbodyptr,
			      rtype*4,&posptr,NDIM)) {
	      fprintf(stderr,"Snap error ### No Positions array\n");
	      exit(1);					
	    }
	    else { /* < get_data_pos */
	      if (SP_io) {
		if (!F_dim) { /* Fortran array : Tab(3,MAXBODY) */
		  for(i=0; i<nBodySelected; i++) {
		    memcpy((char *)(pos_f+i*3*jump),
			   (char *)(posptr+SelectedPart[i]*3*jump),
			   jump3);
		  }						  
		}
		else {		 /* Fortran array : Tab(MAXBODY,3) */
		  for(i=0; i<nBodySelected; i++) {
		    for(k=0; k<3; k++) {
		      memcpy((char *) (pos_f+k*(*size_array)*jump+i*jump),
			     (char *) (posptr+SelectedPart[i]*3*jump+k*jump),
			     jump);
		    }
		  }
		}
	      }
	      else {
		if (!F_dim) { /* Fortran array : Tab(3,MAXBODY) */
		  for(i=0; i<*nbodyptr; i++) {
		    memcpy((char *)(pos_f+i*3*jump),
			   (char *)(posptr+i*3*jump),
			   jump3);
		  }						  
		}
		else {		 /* Fortran array : Tab(MAXBODY,3) */
		  for(i=0; i<*nbodyptr; i++) {
		    for(k=0; k<3; k++) {
		      memcpy((char *) (pos_f+k*(*size_array)*jump+i*jump),
			     (char *) (posptr+i*3*jump+k*jump),
			     jump);
		    }
		  }
		}
	      }
	      free((char *) posptr);
	    } /* > get_data_pos */
	  }
	  if (V_io) {
	    dprintf(1,"Getting Velocities....\n");
	    if (!get_data_vel(instr[no_io],OutType,*nbodyptr,
			      rtype*4,&velptr,NDIM)) {
	      fprintf(stderr,"Snap error ### No Velocities array\n");
	      exit(1);					
	    }
	    else { /* < get_data_vel */
	      if (SP_io) {
		if (!F_dim) { /* Fortran array : Tab(3,MAXBODY) */
		  for(i=0; i<nBodySelected; i++) {
		    memcpy((char *)(vel_f+i*3*jump),
			   (char *)(velptr+SelectedPart[i]*3*jump),
			   jump3);
		  }						  
		}
		else {		 /* Fortran array : Tab(MAXBODY,3) */
		  for(i=0; i<nBodySelected; i++) {
		    for(k=0; k<3; k++) {
		      memcpy((char *) (vel_f+k*(*size_array)*jump+i*jump),
			     (char *) (velptr+SelectedPart[i]*3*jump+k*jump),
			     jump);
		    }
		  }
		}
	      }
	      else {
		if (!F_dim) { /* Fortran array : Tab(3,MAXBODY) */
		  for(i=0; i<*nbodyptr; i++) {
		    memcpy((char *)(vel_f+i*3*jump),
			   (char *)(velptr+i*3*jump),
			   jump3);
		  }						  
		}
		else {		 /* Fortran array : Tab(MAXBODY,3) */
		  for(i=0; i<*nbodyptr; i++) {
		    for(k=0; k<3; k++) {
		      memcpy((char *) (vel_f+k*(*size_array)*jump+i*jump),
			     (char *) (velptr+i*3*jump+k*jump),
			     jump);
		    }
		  }
		}
	      }
	      free((char *) velptr);
	    } /* > get_data_vel */
	  }
	} /* else get_data_pos || get_data_vel */
      } /* if (X_io || V_io  ) { */
      /* get potentials */
      if (P_io) {
	dprintf(1,"Getting Potentials....\n");
	if (!get_data_pot(instr[no_io],OutType,*nbodyptr,
			  rtype*4,&potptr)) {
	  fprintf(stderr,"Snap error ### No Potential\n");
	  exit(1);
	} 	
	else {
	  /*for (i=0; i<*nbodyptr; i++)*/
	  if (SP_io)
	    for (i=0; i<nBodySelected; i++) {
	      memcpy((char*)(pot_f+jump*i),
		     (char*)(potptr+(jump*SelectedPart[i])),
		     jump);
	    }
	  else
	    memcpy((char *) (pot_f),
		   (char *) (potptr),
		   *nbodyptr * jump);		
	  free(potptr);
	}
      } 				 
      /* get accelerations */
      if (A_io) {
	dprintf(1,"Getting Accelerations....\n");
	if (!get_data_acc(instr[no_io],OutType,*nbodyptr,
			  rtype*4,&accptr,NDIM)) {
	  fprintf(stderr,"Snap error ### No Acceleration array\n");
	  exit(1);					
	}
	else {
	  if (SP_io) {
	    if (!F_dim) { /* Fortran array : Tab(3,MAXBODY) */
	      for(i=0; i<nBodySelected; i++) {
		memcpy((char *)(acc_f+i*3*jump),
		       (char *)(accptr+SelectedPart[i]*3*jump),
		       jump3);
	      }						  
	    }
	    else {		 /* Fortran array : Tab(MAXBODY,3) */
	      for(i=0; i<nBodySelected; i++) {
		for(k=0; k<3; k++) {
		  memcpy((char *) (acc_f+k*(*size_array)*jump+i*jump),
			 (char *) (accptr+SelectedPart[i]*3*jump+k*jump),
			 jump);
		}
	      }
	    }
	  }
	  else {
	    if (!F_dim) { /* Fortran array : Tab(3,MAXBODY) */
	      for(i=0; i<*nbodyptr; i++) {
		memcpy((char *)(acc_f+i*3*jump),
		       (char *)(accptr+i*3*jump),
		       jump3);
	      }						  
	    }
	    else {		 /* Fortran array : Tab(MAXBODY,3) */
	      for(i=0; i<*nbodyptr; i++) {
		for(k=0; k<3; k++) {
		  memcpy((char *) (acc_f+k*(*size_array)*jump+i*jump),
			 (char *) (accptr+i*3*jump+k*jump),
			 jump);
		}
	      }
	    }
	  }
	  free((char *) accptr);
	}
      } /* if (A_io) { */
      /* get Aux data */
      if (AUX_io) {
	dprintf(1,"Getting Aux....\n");
	if (!get_data_aux(instr[no_io],OutType,*nbodyptr,
			   rtype*4,&auxptr)) {
	  fprintf(stderr,"Snap error ### No AuxTag\n");
	  exit(1);
	} 	
	else {
	  if (SP_io)
	    for (i=0; i<nBodySelected; i++) {
	      memcpy((char*)(aux_f+i_jump*i),
		     (char*)(auxptr+(i_jump*SelectedPart[i])),
		     i_jump);
	    }
	  else
	    memcpy((char *) (aux_f),
		   (char *) (auxptr),
		   *nbodyptr * i_jump);
						
	  free(auxptr);
	}  
      } /* if (AUX_io) { */

      /* get Keys data */
      if (K_io) {
	dprintf(1,"Getting Keys....\n");
	if (!get_data_keys(instr[no_io],IntType,*nbodyptr,
			   sizeof(int),&keysptr)) {
	  fprintf(stderr,"Snap error ### No KeyTag\n");
	  exit(1);
	} 	
	else {
	  if (SP_io)
	    for (i=0; i<nBodySelected; i++) {
	      memcpy((char*)(keys_f+i_jump*i),
		     (char*)(keysptr+(i_jump*SelectedPart[i])),
		     i_jump);
	    }
	  else
	    memcpy((char *) (keys_f),
		   (char *) (keysptr),
		   *nbodyptr * i_jump);
						
	  free(keysptr);
	}  
      } /* if (K_io) { */
      /* get Dens data */
      if (D_io) {
	dprintf(1,"Getting Density....\n");
	if (!get_data_dens(instr[no_io],OutType,*nbodyptr,
			   rtype*4,&densptr)) {
	  fprintf(stderr,"Snap error ### No DensityTag\n");
	  exit(1);
	} 	
	else {
	  if (SP_io)
	    for (i=0; i<nBodySelected; i++) {
	      memcpy((char*)(dens_f+i_jump*i),
		     (char*)(densptr+(i_jump*SelectedPart[i])),
		     i_jump);
	    }
	  else
	    memcpy((char *) (dens_f),
		   (char *) (densptr),
		   *nbodyptr * i_jump);
						
	  free(densptr);
	}  
      } /* if (D_io) { */
      /* get Eps data */
      if (EPS_io) {
	dprintf(1,"Getting Eps....\n");
	if (!get_data_eps(instr[no_io],OutType,*nbodyptr,
			   rtype*4,&epsptr)) {
	  fprintf(stderr,"Snap error ### No EpsTag\n");
	  exit(1);
	} 	
	else {
	  if (SP_io)
	    for (i=0; i<nBodySelected; i++) {
	      memcpy((char*)(eps_f+i_jump*i),
		     (char*)(epsptr+(i_jump*SelectedPart[i])),
		     i_jump);
	    }
	  else
	    memcpy((char *) (eps_f),
		   (char *) (epsptr),
		   *nbodyptr * i_jump);
						
	  free(epsptr);
	}  
      } /* if (EPS_io) { */
      get_tes(instr[no_io], ParticlesTag);
							
      /* get out of the loop */
      get_tes(instr[no_io], SnapShotTag);
      break;
    }
    else {
      fprintf(stderr,"Snap error ### no Particles\n");
      exit(1);
    } 
    
    /* get_tes(instr[no_io], SnapShotTag); */
  }
  /* free allocated memory */
  if (ST_io) {
    free((char *) select_time);
  }
  if (SP_io) {
    free((int *) SelectedPart);
    free((char*) SelString);
    *nbody_f=nBodySelected;
  }
  /* check out dynamic snapshots */
  if (maxbodies[CURRENT_IO] < *nbodyptr) {
    maxbodies[CURRENT_IO] = *nbodyptr;
  }
  /* garbage collecting */
  free((int *) nbodyptr); 
  free((char *) timeptr);

  dprintf(1,"End of [get_data_select_f]....\n");
  return 1;
}
/* ----------------------------------------------------------------
|  End of io_get_put_f.c                                           
+---------------------------------------------------------------- */
