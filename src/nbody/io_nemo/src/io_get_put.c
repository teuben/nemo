/* -------------------------------------------------------------- *\
|* $Id$
|* 
|* Get/Put selected Data from io_nemo
|*
\* -------------------------------------------------------------- */

/* -------------------------------------------------------------- *\
|* Nemo include files
\* -------------------------------------------------------------- */
#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>		
#include <filestruct.h>

#include <snapshot/snapshot.h>	
#include <snapshot/body.h>
#include <history.h>

/* -------------------------------------------------------------- *\
|* Local include files
\* -------------------------------------------------------------- */
#include "io_get_put.h"
#include "check_file.h"
#include "io_nemo_tools.h"
#include "get_data_wrapper.h"

/* extern variables */
#include "history_data.h"
#include "flags_data.h"
#include "io_nemo_data.h"


#define TIMEFUZZ 0.0000001

#ifdef __cplusplus
extern "C" {
#endif
  bool within();
	
#ifdef __cplusplus
}
#endif
/* -------------------------------------------------------------- *\ 
|* put_data_select :
|* Save the snapshot in NEMO format.
\* -------------------------------------------------------------- */
int put_data_select(char * outfile,
		    int    rtype,
		    char * io_out[],
		    bool * save_one,
		    FILE * outstr[],
		    int    MAXIO)
{
  int coordsys =    CSCode(Cartesian, NDIM, 2),
    i,j,k,no_io;
  char * phasep;    /* Phase space coordinates */
  /* char * accptr; */   /* Acceleration array */
	
  int jump = rtype*4,jump3=3*jump,jump6=6*jump;
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
    phasep = (char *) malloc(sizeof(char)*2*NDIM*(*nbody)*jump);

    if (phasep == NULL) {
      fprintf(stderr,"Memory error ## [put_data_select]\n");
      fprintf(stderr,"Impossible to allocate PhaseSpace Coordinates\n");
    }
		
    /* store positions and velocities in phasep pointer */
    for(i=0; i<*nbody; i++) {
      if (X_io)
	memcpy((char *) (phasep+i*jump6),
	       (char *) (pos+i*jump3),jump3);
			
      if (V_io)
	memcpy((char *) (phasep+i*jump6+jump3),
	       (char *) (vel+i*jump3),
	       jump3);
    }  
  }
#endif
  /* set up the history */
  if (!set_history[no_io]) {
    /* set TRUE to history */
    set_history[no_io]=TRUE;
		
    /* raz history */
    reset_history();
		
    /* build command line's history */
    app_history(history_prog); 
		
    if (H_io) { /* request for history file name */
      /* get history from input file */
      get_history_input_file(hist_file);
    }
    /* put the history */
    put_history(outstr[no_io]);
  }
	
  /* Save the snapshot  */ 
  put_set(outstr[no_io], SnapShotTag);	
	
  put_set(outstr[no_io], ParametersTag);
  if (T_io)
    put_data(outstr[no_io], TimeTag, OutType, timu, 0);
  put_data(outstr[no_io], NobjTag, IntType, nbody, 0);
  put_tes(outstr[no_io], ParametersTag);
	
  put_set(outstr[no_io], ParticlesTag); /*   start particle output  */
  put_data(outstr[no_io], CoordSystemTag,IntType,&coordsys,0);

  if (M_io)      /* Masses */
    put_data(outstr[no_io],MassTag,OutType,mass,*nbody,0);

  if (XV_io) {   /* PhaseSpace */
    put_data(outstr[no_io],PhaseSpaceTag,OutType,phase,*nbody,2,NDIM,0); 
    free(phasep); 
  }

  if (X_io)     /* Positions */
    put_data(outstr[no_io],PosTag,OutType,pos,*nbody,NDIM,0);

  if (V_io)     /* Velocities */
    put_data(outstr[no_io],VelTag,OutType,vel,*nbody,NDIM,0);

  if (P_io)     /* Potentials */
    put_data(outstr[no_io], PotentialTag,OutType,pot,*nbody,0);

  if (A_io)     /* Accelerations */
    put_data(outstr[no_io], AccelerationTag,OutType,acc,*nbody,NDIM, 0);

  if (K_io)     /* Keys */
    put_data(outstr[no_io],KeyTag,IntType,keys,*nbody,0);

  put_tes(outstr[no_io], ParticlesTag);
	
  put_tes(outstr[no_io], SnapShotTag);
	
  /* flush IO */
  fflush(outstr[no_io]);
	
  save_one[no_io] = TRUE;
	
  return 1;
	
}
/* -------------------------------------------------------------- *\ 
|* get_data_select :
|* Read the snapshot in NEMO format.
\* -------------------------------------------------------------- */
int get_data_select(char * infile,
		    int    rtype,/* real type : 1=real*4 | 2=real*8 */
		    char * io_in[],
		    bool   read_one[],
		    FILE * instr[],
		    int    MAXIO) 
{
  char * phaseptr = NULL;
  int  * nbodyptr = NULL;
  int i,j,k, no_io,i_jump,jump=rtype*4,jump3=3*jump,jump6=6*jump;
    
  string headline, * histo;
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
      fprintf(stderr,"WARNING!! end of snapshot reached\n");
      read_one[no_io] = FALSE;
      free(io_in[no_io]);
      strclose(instr[no_io]);
      return 0;
    } /* !get_tag_ok(instr[no_io], SnapShotTag)... */
    else {
      read_one[no_io] = TRUE; /* File has been read at least one time */
    }		
    get_set(instr[no_io], SnapShotTag);
		
    get_set(instr[no_io], ParametersTag);
    if (T_io) {  /* time step */
      get_data_time(instr[no_io], OutType, rtype*4, &timu);
			
    }

    get_data_nbody(instr[no_io], IntType, sizeof(int),/*(int*)*/ &nbodyptr);

    if (SP_io) {
      if ( strcmp(SelectionString123,"all")) {
	SelectedPart=(int*) allocate(*nbodyptr*sizeof(int));
	/*create index list */
	nBodySelected=nemoinpi(SelectionString123,
			       SelectedPart,*nbodyptr); 
	if (nBodySelected < 0) {
	  fprintf(stderr,"Failed to select particles's range <%s>"
		  " *nemoinpi* function return code = [%d],"
		  " aborted.....\n",SelectionString123,nBodySelected);
	  exit(1);
	}
      }
      else { /* all particles are selected, so we *DESACTVATE* SP_io */
	SP_io = 0;
      }
    }
    if (N_io) { /* nbody */
      nbody = (int *) allocate_pointer((char *)nbody,4);
      *nbody = *nbodyptr;
    }
    if (ST_io) { /* selected time */
				/* get the real value of the selected time */
      real_time= char2double(timu,rtype);
			
      if (!streq(selt,"all") && !within(real_time,selt,TIMEFUZZ)) {
	fprintf(stderr,"Info : skipping time step [%.4f]\n",real_time);
	get_tes(instr[no_io], ParametersTag);
	get_tes(instr[no_io], SnapShotTag);
	continue; /* go to the next time step */
      }
    }
      
    get_tes(instr[no_io], ParametersTag);

    if (get_tag_ok(instr[no_io], ParticlesTag)) {
      get_set(instr[no_io], ParticlesTag);
			
      /* get masses */
      if (M_io)
	if (!get_data_mass(instr[no_io],OutType,*nbodyptr,rtype*4,&mass)) {
	  fprintf(stderr,"Snap error ### No Mass\n");
	  exit(1);
	}
	else 
	  if (SP_io) {
	    for (i=0;i<nBodySelected;i++) {
	      memcpy ((char *) (mass+i*jump),
		      (char *) (mass+SelectedPart[i]*jump),
		      jump);
	    }
	  }

      /* get positions and velocities */
      if (X_io || V_io || XV_io) {

	/* test PHASE SPACE Coordinates */
	if ( get_data_phase(instr[no_io],OutType,*nbodyptr, rtype*4,&phaseptr,NDIM)) {
	  
	  if (X_io)
	    pos = (char *) allocate_pointer(pos,*nbody*3*4*rtype);
	  if (V_io)
	    vel = (char *) allocate_pointer(vel,*nbody*3*4*rtype);
					
	  if (SP_io) {
	    if (X_io)
	      for (i=0;i<nBodySelected;i++) {
		memcpy((char *) (pos+i*jump3),
		       (char *) (phaseptr+SelectedPart[i]*jump6),
		       jump3);
	      }
	    if (V_io)
	      for (i=0;i<nBodySelected;i++) {
		memcpy((char *) (vel+i*jump3),
		       (char *) (phaseptr+SelectedPart[i]*jump6+jump3),
		       jump3);
	      }
	    if (XV_io)
	      for (i=0;i<nBodySelected;i++) {
		memcpy((char *) (phase+i*jump6),
		       (char *) (phaseptr+SelectedPart[i]*jump6),
		       jump6);
	      }
	  }
	  else {  /* ! SP_io */
	    
	    for(i=0; i<*nbodyptr; i++) {
		if (X_io)
		  memcpy((char *) (pos+i*jump3),
			 (char *) (phaseptr+i*jump6),
			 jump3);
		if (V_io)
		  memcpy((char *) (vel+i*jump3),
			 (char *) (phaseptr+i*jump6+jump3),
			 jump3);		
	    } /* for */                         
	    if (XV_io)
	      memcpy((char *) phase,
		     (char *) phaseptr,
		     (*nbodyptr)*jump6);
	  }
	  free(phaseptr);

	} /* if (get_data_phase */
	
	else { /* PosTag and VelTag */
	  
	  if (X_io) {

	    /* test PosTag */
	    if ( get_data_pos(instr[no_io],OutType,*nbodyptr, rtype*4,&pos,NDIM)) {
	      if (SP_io) {
		for (i=0;i<nBodySelected;i++) {
		  memcpy ((char *) (pos+i*jump3),
			  (char *) (pos+SelectedPart[i]*jump3),
			  jump3);
		}
	      }
	    }
	    else {
	      fprintf(stderr,"Snap error ### No Positions\n");
	      exit(1);
	    }
	  } /* if (X_io) */

	  if (V_io) {
	    /* test VelTag */
	    if ( get_data_vel(instr[no_io],OutType,*nbodyptr, rtype*4,&vel,NDIM)) {
	      if (SP_io) {
		for (i=0;i<nBodySelected;i++) {
		  memcpy ((char *) (vel+i*jump3),
			  (char *) (vel+SelectedPart[i]*jump3),
			  jump3);
		}
	      }
	    }
	    else {
	      fprintf(stderr,"Snap error ### No Velocities\n");
	      exit(1);
	    }
	  } /* if (V_io)             */
	} /* else PosTag and VelTag */
      } /* if  (X_io || V_io)      */
    
      /* get potential array */
      if (P_io)
	if (!get_data_pot(instr[no_io],OutType,*nbodyptr, rtype*4,&pot)) {
	  fprintf(stderr,"Snap error ### No Potential\n");
	  exit(1);
	}  
	else 
	  if (SP_io) {
	    for (i=0;i<nBodySelected;i++) {
	      memcpy ((char *) (pot+i*jump),
		      (char *) (pot+SelectedPart[i]*jump),
		      jump);
	    }
	  }
      /* get acceleration array */
      if (A_io)
	if (!get_data_acc(instr[no_io],OutType,*nbodyptr, rtype*4,&acc,NDIM)) {
	  fprintf(stderr,"Snap error ### No Acceleration\n");
	  exit(1);                    
	}
	else 
	  if (SP_io) {
	    for (i=0;i<nBodySelected;i++) {
	      memcpy((char *) (acc+i*jump3),
		     (char *) (acc+SelectedPart[i]*jump3),
		     jump3);
	    }
	  }
      /* get keys array */
      if (K_io)
	if (!get_data_keys(instr[no_io],IntType,*nbodyptr, rtype*4,&keys)) {
	  fprintf(stderr,"Snap error ### No Keys\n");
	  exit(1);
	}  
	else 
	  if (SP_io) {
	    for (i=0;i<nBodySelected;i++) {
	      memcpy ((char *) (keys+i*i_jump),
		      (char *) (keys+SelectedPart[i]*i_jump),
		      i_jump);
	    }
	  }
      get_tes(instr[no_io], ParticlesTag);
    } /*  get_tag_ok(instr[no_io], ParticlesTag */
    else {
      fprintf(stderr,"Snap error ### no Particles\n");
      exit(1);
    } 
    
    /* get out of the loop */
    get_tes(instr[no_io], SnapShotTag);
    break;								
  } /* for(;;) */

  /* get_tes(instr[no_io], SnapShotTag); */
  free(nbodyptr);
  
  if (SP_io) {
    free(SelectedPart);
    free(SelectionString123);
    *nbody=nBodySelected;
  }
  return 1;

}
/* -------------------------------------------------------------- *\
|* End of io_get_put.c
\* -------------------------------------------------------------- */
