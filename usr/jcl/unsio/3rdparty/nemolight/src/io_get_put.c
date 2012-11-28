/* =================================================================
|  Copyright Jean-Charles LAMBERT - 2008                            
|  e-mail:   Jean-Charles.Lambert@oamp.fr                           
|  address:  Dynamique des galaxies                                 
|            Laboratoire d'Astrophysique de Marseille               
|            2, place Le Verrier                                    
|            13248 Marseille Cedex 4, France                        
|            CNRS U.M.R 6110                                        
| ==================================================================
|* Get/Put selected Data from io_nemo                               
| ==================================================================
| 03-Mar-05 : a) critical memory bug corrected which could corrupt  
|                memory in case of using io_nemo with several       
|                snapshots.                                         
|             b) nemo keybits added.                                
+----------------------------------------------------------------- */

/* -----------------------------------------------------------------
|  Nemo include files                                               
+----------------------------------------------------------------- */
#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>		
#include <filestruct.h>

#include <snapshot/snapshot.h>	
#include <snapshot/body.h>
#include <history.h>

/* ----------------------------------------------------------------
|  Local include files                                             
+---------------------------------------------------------------- */
#include "io_get_put.h"
#include "check_file.h"
#include "io_nemo_tools.h"
#include "get_data_wrapper.h"

/* extern variables */
#include "history_data.h"
#include "flags_data.h"
#include "io_nemo_data.h"
#include "parameters.h"


#define TIMEFUZZ 0.0000001
extern int maxbodies[];
int CURRENT_IO;
/* -----------------------------------------------------------------
|  put_data_select :                                                
|  Save the snapshot in NEMO format.                                
+----------------------------------------------------------------- */
int put_data_select(char * outfile,
		    int    rtype,
		    char * io_out[],
		    bool * save_one,
		    FILE * outstr[],
		    int    MAXIO,
		    t_ion_data * ion)
{
  int coordsys =    CSCode(Cartesian, NDIM, 2),
    no_io;
  /*char * phasep;*/    /* Phase space coordinates */
  /* char * accptr; */   /* Acceleration array */
	
  /*int jump = rtype*4;*/
  char *  OutType;
	
  if (rtype==1)
    OutType=FloatType;
  else
    OutType=DoubleType;
	
	
  /* Open the file for writing */
  if ((no_io = get_old_file(outfile,io_out,save_one,outstr,MAXIO)) < 0)
    no_io = get_new_file(outfile,io_out,save_one,outstr,"w",MAXIO);
	
  /* print out what it is doing */
  if (I_io)
    chk_parameters(FALSE,0, rtype);
	
  /* memory allocation for positions and velocities */
#if 0
  if (XV_io) {
    phasep = (char *) malloc(sizeof(char)*2*NDIM*(*ion->nbody)*jump);

    if (phasep == NULL) {
      fprintf(stderr,"Memory error ## [put_data_select]\n");
      fprintf(stderr,"Impossible to allocate PhaseSpace Coordinates\n");
    }
		
    /* store positions and velocities in phasep pointer */
    for(i=0; i<*ion->nbody; i++) {
      if (X_io)
	memcpy((char *) (phasep+i*jump6),
	       (char *) (ion->pos+i*jump3),jump3);
			
      if (V_io)
	memcpy((char *) (phasep+i*jump6+jump3),
	       (char *) (ion->vel+i*jump3),
	       jump3);
    }  
  }
#endif

  /* set up the history */
  if (!set_history[no_io]) {
    /* set TRUE to history */
    set_history[no_io]=TRUE;
		
    /* raz history */
    if (!( H_io  && ! strcmp(hist_file,"-"))) { // do not reset history if hist_file = "-"
      reset_history();
      /* build command line's history */
      app_history(history_prog); 
    }
		
    if (H_io && strcmp(hist_file,"-")) { /* request for history file name */
      /* get history from input file */
      get_history_input_file(hist_file);
    }
    /* put the history */
    put_history(outstr[no_io]);
  }

  /* Save the snapshot  */ 
  put_set(outstr[no_io], SnapShotTag);	
	
  put_set(outstr[no_io], ParametersTag);
  if (T_io) {
    if ( (B_io && ( *ion->bits & TimeBit)) || !B_io) {
      put_data(outstr[no_io], TimeTag, OutType, ion->timu, 0);
    } 
    else {
      dprintf(1,"WARNING ### TimeBit control does not exist.\n");
    }
  }
  put_data(outstr[no_io], NobjTag, IntType, ion->nbody, 0);
  put_tes(outstr[no_io], ParametersTag);
	
  put_set(outstr[no_io], ParticlesTag); /*   start particle output  */
  put_data(outstr[no_io], CoordSystemTag,IntType,&coordsys,0);

  if (M_io) {      /* Masses */
    if ( (B_io && ( *ion->bits & MassBit)) || !B_io) {
      put_data(outstr[no_io],MassTag,OutType,ion->mass,*ion->nbody,0);
    } 
    else {
      dprintf(1,"WARNING ### MassBit control does not exist.\n");
    }
  }

  if (XV_io) {   /* PhaseSpace */
    if ( (B_io && ( *ion->bits & PhaseSpaceBit)) || !B_io) {
      put_data(outstr[no_io],PhaseSpaceTag,OutType,ion->phase,*ion->nbody,2,NDIM,0); 
      /* free(phasep); */
    }
    else {
      dprintf(1,"WARNING ### PhaseSpaceBit control does not exist.\n");
    }
  }

  if (X_io) {     /* Positions */
    if ( (B_io && ( *ion->bits & PosBit)) || !B_io) {
      put_data(outstr[no_io],PosTag,OutType,ion->pos,*ion->nbody,NDIM,0);
    }
    else {
      dprintf(1,"WARNING ### PosBit control does not exist.\n");
    }
  }

  if (V_io) {     /* Velocities */
    if ( (B_io && ( *ion->bits & VelBit)) || !B_io) {
      put_data(outstr[no_io],VelTag,OutType,ion->vel,*ion->nbody,NDIM,0);
    }
    else {
      dprintf(1,"WARNING ### VelBit control does not exist.\n");
    }
  }

  if (P_io) {     /* Potentials */
    if ( (B_io && ( *ion->bits & PotentialBit)) || !B_io) {
      put_data(outstr[no_io], PotentialTag,OutType,ion->pot,*ion->nbody,0);
    }
    else {
      dprintf(1,"WARNING ### PotentialBit control does not exist.\n");
    }
  }

  if (A_io) {     /* Accelerations */
    if ( (B_io && ( *ion->bits & AccelerationBit)) || !B_io) {
      put_data(outstr[no_io], AccelerationTag,OutType,ion->acc,*ion->nbody,NDIM, 0);
    }
    else {
      dprintf(1,"WARNING ### AccelerationBit control does not exist.\n");
    }
  }

  if (AUX_io) {      /* Aux */
    if ( (B_io && ( *ion->bits & AuxBit)) || !B_io) {
      put_data(outstr[no_io],AuxTag,OutType,ion->aux,*ion->nbody,0);
    } 
    else {
      dprintf(1,"WARNING ### AuxBit control does not exist.\n");
    }
  }

  if (K_io) {     /* Keys */
    if ( (B_io && ( *ion->bits & KeyBit)) || !B_io) {
      put_data(outstr[no_io],KeyTag,IntType,ion->keys,*ion->nbody,0);
    }
    else {
      dprintf(1,"WARNING ### KeyBit control does not exist.\n");
    }
  }

  if (D_io) {      /* Density */
    if ( (B_io && ( *ion->bits & DensBit)) || !B_io) {
      put_data(outstr[no_io],DensityTag,OutType,ion->dens,*ion->nbody,0);
    } 
    else {
      dprintf(1,"WARNING ### DensBit control does not exist.\n");
    }
  }

  if (EPS_io) {   /* EPS */
    if ( (B_io && ( *ion->bits & EpsBit)) || !B_io) {
      put_data(outstr[no_io],EpsTag,OutType,ion->eps,*ion->nbody,0);
    }
    else {
      dprintf(1,"WARNING ### EpsBit control does not exist.\n");
    }
  }

  put_tes(outstr[no_io], ParticlesTag);
	
  put_tes(outstr[no_io], SnapShotTag);
	
  /* flush IO */
  fflush(outstr[no_io]);
	
  save_one[no_io] = TRUE;
	
  return 1;
	
}
/* -----------------------------------------------------------------
|  get_data_select :                                                
|  Read the snapshot in NEMO format.                                
+----------------------------------------------------------------- */
int get_data_select(char * infile,
		    int    rtype,/* real type : 1=real*4 | 2=real*8 */
		    char * io_in[],
		    bool   read_one[],
		    FILE * instr[],
		    int    MAXIO,
		    t_ion_data * ion) 
{
  int keybits=0;   /* store control bits     */
  int status=1;    /* return function status */
  char * phaseptr = NULL;
  int  * nbodyptr = NULL;
  int i, no_io,i_jump,jump=rtype*4,jump3=3*jump,jump6=6*jump;
    
  string headline;
  char *  OutType;
	
  double real_time; /*, char2double();*/
  /*bool within();*/

  int * SelectedPart=NULL; /* added, table of indexes */
  int nBodySelected=0;

  if (rtype==1)
    OutType=FloatType;
  else
    OutType=DoubleType;
	
  i_jump = sizeof(int);
	
  /* Open the file for reading */
  if ((no_io = get_old_file(infile,io_in,read_one,instr,MAXIO)) < 0) {
    no_io = get_new_file(infile,io_in,read_one,instr,"r",MAXIO);
  }
  CURRENT_IO = no_io;
  
  /* print out what is doing */
  if (I_io)
    chk_parameters(TRUE,0, rtype);
	
  /* read the snapshot */
#if 0
  reset_history();
  fprintf(stderr,"nhist : %d\n",get_history(instr[no_io]));    
  histo = (char **) ask_history();
  for (i=0; i<get_history(instr[no_io]); i++) {
    fprintf(stderr,"histo : <%s>\n",histo[i]);
  }
#endif
  reset_history();
  for (;;) {

    get_history(instr[no_io]);

    while (get_tag_ok(instr[no_io],HeadlineTag))
      headline = get_string(instr[no_io],HeadlineTag);

    if (!get_tag_ok(instr[no_io], SnapShotTag)) {
      if (!read_one[no_io]) {
	fprintf(stderr,"SnapshotTag error ## [get_data_select]\n");
	fprintf(stderr,"%s is not a NEMO SNAPSHOT\n",io_in[no_io]);
	exit(1);
      }
      /* end of snapshot reached */
      dprintf(1,"WARNING!! end of snapshot reached.\n");

      return 0;
    } /* !get_tag_ok(instr[no_io], SnapShotTag)... */
    else {
      read_one[no_io] = TRUE; /* File has been read at least one time */
    }		
    get_set(instr[no_io], SnapShotTag);
		
    get_set(instr[no_io], ParametersTag);
    if (T_io) {  /* time step */
      if (!get_data_time(instr[no_io], OutType, rtype*4, &ion->timu)) {
	dprintf(1,"### Snapshot WARNING ### No Time\n");
	status = -1;
      }
      else {
	keybits |= TimeBit;   /* got Time */
      }
			
    }

    get_data_nbody(instr[no_io], IntType, sizeof(int),/*(int*)*/ &nbodyptr);

    if (SP_io) {
      if ( strcmp(ion->SelectionString123,"all")) {
	SelectedPart=(int*) allocate(*nbodyptr*sizeof(int));
	/*create index list */
	nBodySelected=nemoinpi(ion->SelectionString123,
			       SelectedPart,*nbodyptr); 
	if (nBodySelected < 0) {
	  fprintf(stderr,"Failed to select particles's range <%s>"
		  " *nemoinpi* function return code = [%d],"
		  " aborted.....\n",ion->SelectionString123,nBodySelected);
	  exit(1);
	}
      }
      else { /* all particles are selected, so we *DESACTVATE* SP_io */
	SP_io = 0;
      }
    }
    if (N_io) { /* nbody */
      ion->nbody = (int *) allocate_pointer((char *)ion->nbody,4);
      *ion->nbody = *nbodyptr;
    }
    if (ST_io && (keybits & TimeBit) ) { /* selected time */
				/* get the real value of the selected time */
      real_time= char2double(ion->timu,rtype);
			
      if (!streq(ion->selt,"all") && !within(real_time,ion->selt,TIMEFUZZ)) {
	dprintf(1,"Info : skipping time step [%.4f]\n",real_time);
	get_tes(instr[no_io], ParametersTag);
	get_tes(instr[no_io], SnapShotTag);
	if (SP_io) {  /* stuff on Selected Particles */
	  free((int *) SelectedPart);
	}
	continue; /* go to the next time step */
      }
    }
      
    get_tes(instr[no_io], ParametersTag);

    if (get_tag_ok(instr[no_io], ParticlesTag)) {
      get_set(instr[no_io], ParticlesTag);
			
      /* get masses */
      if (M_io) {
	if (!get_data_mass(instr[no_io],OutType,*nbodyptr,rtype*4,&ion->mass)) {
	  dprintf(1,"### Snapshot WARNING ### No Mass\n");
	  status = -1;
	  /*exit(1);*/
	}
	else {
	  keybits |= MassBit;   /* got Mass */
	  if (SP_io) {
	    for (i=0;i<nBodySelected;i++) {
	      memcpy ((char *) (ion->mass+i*jump),
		      (char *) (ion->mass+SelectedPart[i]*jump),
		      jump);
	    }
	  }
	}
      }
      /* get positions and velocities */
      if (X_io || V_io || XV_io) {

	/* test PHASE SPACE Coordinates */
	if ( get_data_phase(instr[no_io],OutType,*nbodyptr, rtype*4,&phaseptr,NDIM)) {
	  keybits |= PhaseSpaceBit;   /* got PhaseSpace */
	  if (X_io) {
	    keybits |= PosBit;   /* got Pos */
	    if (maxbodies[CURRENT_IO] < *nbodyptr) {
	      if (ion->pos) {
		free ((char *) (ion->pos));
		ion->pos = NULL;
	      }
	      
	    }
	    ion->pos = (char *) allocate_pointer(ion->pos,*ion->nbody*3*4*rtype);
	  }
	  if (V_io) {
	    keybits |= VelBit;   /* got Vel */
	    if (maxbodies[CURRENT_IO] < *nbodyptr) {
	      if (ion->vel) {
		free ((char *) (ion->vel));
		ion->vel = NULL;
	      }
	    }
	    ion->vel = (char *) allocate_pointer(ion->vel,*ion->nbody*3*4*rtype);
	  }
					
	  if (SP_io) {
	    if (X_io)
	      for (i=0;i<nBodySelected;i++) {
		memcpy((char *) (ion->pos+i*jump3),
		       (char *) (phaseptr+SelectedPart[i]*jump6),
		       jump3);
	      }
	    if (V_io)
	      for (i=0;i<nBodySelected;i++) {
		memcpy((char *) (ion->vel+i*jump3),
		       (char *) (phaseptr+SelectedPart[i]*jump6+jump3),
		       jump3);
	      }
	    if (XV_io)
	      for (i=0;i<nBodySelected;i++) {
		memcpy((char *) (ion->phase+i*jump6),
		       (char *) (phaseptr+SelectedPart[i]*jump6),
		       jump6);
	      }
	  }
	  else {  /* ! SP_io */
	    
	    for(i=0; i<*nbodyptr; i++) {
		if (X_io)
		  memcpy((char *) (ion->pos+i*jump3),
			 (char *) (phaseptr+i*jump6),
			 jump3);
		if (V_io)
		  memcpy((char *) (ion->vel+i*jump3),
			 (char *) (phaseptr+i*jump6+jump3),
			 jump3);		
	    } /* for */                         
	    if (XV_io)
	      memcpy((char *) ion->phase,
		     (char *) phaseptr,
		     (*nbodyptr)*jump6);
	  }
	  free(phaseptr);

	} /* if (get_data_phase */
	else { /* PosTag and VelTag */
	  
	  if (X_io) {

	    /* test PosTag */
	    if ( get_data_pos(instr[no_io],OutType,*nbodyptr, rtype*4,&(ion->pos),NDIM)) {
	      keybits |= PosBit;   /* got Pos */
	      if (SP_io) {
		for (i=0;i<nBodySelected;i++) {
		  memcpy ((char *) (ion->pos+i*jump3),
			  (char *) (ion->pos+SelectedPart[i]*jump3),
			  jump3);
		}
	      }
	    }
	    else {
	      dprintf(1,"### Snapshot WARNING ### No Positions\n");
	      status = -1;
	      /*exit(1);*/
	    }
	  } /* if (X_io) */

	  if (V_io) {
	    /* test VelTag */
	    if ( get_data_vel(instr[no_io],OutType,*nbodyptr, rtype*4,&ion->vel,NDIM)) {
	      keybits |= VelBit;   /* got Vel */
	      if (SP_io) {
		for (i=0;i<nBodySelected;i++) {
		  memcpy ((char *) (ion->vel+i*jump3),
			  (char *) (ion->vel+SelectedPart[i]*jump3),
			  jump3);
		}
	      }
	    }
	    else {
	      dprintf(1,"### Snapshot WARNING ### No Velocities\n");
	      status = -1;
	      /*exit(1);*/
	    }
	  } /* if (V_io)             */
	} /* else PosTag and VelTag */
      } /* if  (X_io || V_io)      */
    
      /* get potential array */
      if (P_io) {
	if (!get_data_pot(instr[no_io],OutType,*nbodyptr, rtype*4,&ion->pot)) {
	  dprintf(1,"### Snapshot WARNING ### No Potential\n");
	  status = -1;
	  /*exit(1);*/
	}  
	else {
	  keybits |= PotentialBit;   /* got Potential */
	  if (SP_io) {
	    for (i=0;i<nBodySelected;i++) {
	      memcpy ((char *) (ion->pot+i*jump),
		      (char *) (ion->pot+SelectedPart[i]*jump),
		      jump);
	    }
	  }
	}
      }
      /* get acceleration array */
      if (A_io) {
	if (!get_data_acc(instr[no_io],OutType,*nbodyptr, rtype*4,&ion->acc,NDIM)) {
	  dprintf(1,"### Snapshot WARNING ### No Acceleration\n");
	  status = -1;
	  /*exit(1);*/
	}
	else {
	  keybits |= AccelerationBit;   /* got Pos */
	  if (SP_io) {
	    for (i=0;i<nBodySelected;i++) {
	      memcpy((char *) (ion->acc+i*jump3),
		     (char *) (ion->acc+SelectedPart[i]*jump3),
		     jump3);
	    }
	  }
	}
      }
      /* get auxiliary array */
      if (AUX_io) {
	if (!get_data_aux(instr[no_io],OutType,*nbodyptr, rtype*4,&ion->aux)) {
	  dprintf(1,"### Snapshot WARNING ### No Auxiliary\n");
	  status = -1;
	  /*exit(1);*/
	}  
	else {
	  keybits |= AuxBit;   /* got Auxiliary */
	  if (SP_io) {
	    for (i=0;i<nBodySelected;i++) {
	      memcpy ((char *) (ion->aux+i*jump),
		      (char *) (ion->aux+SelectedPart[i]*jump),
		      jump);
	    }
	  }
	}
      }
      /* get keys array */
      if (K_io) {
	if (!get_data_keys(instr[no_io],IntType,*nbodyptr, rtype*4,&ion->keys)) {
	  dprintf(1,"### Snapshot WARNING ### No Keys\n");
	  status = -1;
	  /*exit(1);*/
	}  
	else {
	  keybits |= KeyBit;   /* got Keys */
	  if (SP_io) {
	    for (i=0;i<nBodySelected;i++) {
	      memcpy ((char *) (ion->keys+i*i_jump),
		      (char *) (ion->keys+SelectedPart[i]*i_jump),
		      i_jump);
	    }
	  }
	}
      }
      /* get density array */
      if (D_io) {
	if (!get_data_dens(instr[no_io],OutType,*nbodyptr, rtype*4,&ion->dens)) {
	  dprintf(1,"### Snapshot WARNING ### No Density\n");
	  status = -1;
	  /*exit(1);*/
	}  
	else {
	  keybits |= DensBit;   /* got Density */
	  if (SP_io) {
	    for (i=0;i<nBodySelected;i++) {
	      memcpy ((char *) (ion->dens+i*jump),
		      (char *) (ion->dens+SelectedPart[i]*jump),
		      jump);
	    }
	  }
	}
      }
      /* get epses */
      if (EPS_io) {
	if (!get_data_eps(instr[no_io],OutType,*nbodyptr,rtype*4,&ion->eps)) {
	  dprintf(1,"### Snapshot WARNING ### No Eps\n");
	  status = -1;
	  /*exit(1);*/
	}
	else {
	  keybits |= EpsBit;   /* got Eps */
	  if (SP_io) {
	    for (i=0;i<nBodySelected;i++) {
	      memcpy ((char *) (ion->eps+i*jump),
		      (char *) (ion->eps+SelectedPart[i]*jump),
		      jump);
	    }
	  }
	}
      }
      get_tes(instr[no_io], ParticlesTag);
    } /*  get_tag_ok(instr[no_io], ParticlesTag */
    else {
      dprintf(1,"### Snapshot WARNING ### no ParticlesTag\n");
      status = -2;
    } 
    
    /* get out of the loop */
    get_tes(instr[no_io], SnapShotTag);
    break;								
  } /* for(;;) */

  if (SP_io) {  /* stuff on Selected Particles */
    free(SelectedPart);
    *ion->nbody=nBodySelected;
  }
  
  if (B_io) {   /* stuff on snapshot control bit */
    ion->bits = (int *) allocate_pointer((char *) ion->bits,sizeof(int));
    *ion->bits = keybits;
  }
  /* check out dynamic snapshots */
  if (maxbodies[CURRENT_IO] < *nbodyptr) {
    maxbodies[CURRENT_IO] = *nbodyptr;
  }

  /* get_tes(instr[no_io], SnapShotTag); */
  free(nbodyptr);

  return status;

}
/* ----------------------------------------------------------------
|  End of io_get_put.c                                             
+---------------------------------------------------------------- */
