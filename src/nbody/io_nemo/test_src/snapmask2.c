/* -------------------------------------------------------------- *\
|*
|* snapmask2.c :	JC LAMBERT	 14-Mar-95	V1.0
|*					 15-Nov-95      V1.1
|*                      PJT              14-Jul-01      V2.0 (Bastille day!)
|*                      JCL              09-Jul-02      V2.1
|*
|*
|* Copy selected particles at a selected time quickly using few
|* memory.
|*
|* V1.0   : created
|* V1.1   : minor fixes
|* V2.0   : imported in teuben's nemo CVS tree
|* V2.1   : PosTag VelTag and KeyTag added
|* V2.1a  : note that within() uses real, not double

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
#include <history.h>
#define REALLOC

/* -------------------------------------------------------------- *\
|* variables globales
\* -------------------------------------------------------------- */

string defv[] = {		/* DEFAULT INPUT PARAMETERS */
    "in=???\n			  input  NEMO snapshot ",
    "out=???\n			  output NEMO snapshot",
    "select=all\n		  selected particles",
    "times=all\n		  selected times",
    "step=1\n			  copy only step by step snapshots",
    "VERSION=2.1\n		  09-Jul-02 JCL",
    NULL,
};

string usage = "mask out particles while copying snapshot data";

#define TIMEFUZZ	0.00001
/* -------------------------------------------------------------- *\ 
|* save_parameter :
|* Save ParametersTag
\* -------------------------------------------------------------- */
int save_parameter(outstr,time_d,nbody)
     stream outstr;
     double time_d;
     int  nbody;
{
  real time = time_d;
  put_set(outstr, SnapShotTag);		
  put_set(outstr, ParametersTag);
  put_data(outstr, TimeTag, RealType, &time, 0);
  put_data(outstr, NobjTag, IntType, &nbody, 0);
  put_tes(outstr, ParametersTag);  
  return 1;
}
/* -------------------------------------------------------------- *\ 
|* copy_onedim :
|* Copy selected particles from an one dimensionnal real array
\* -------------------------------------------------------------- */
int copy_onedim(Qsel,oneptr,nbody)
     bool * Qsel[];
     real * oneptr;
     int nbody;
{
  real * op, * dp;
  int  i;

  for (i=0, op=oneptr, dp = oneptr; i<nbody; i++, op++) {
    if (!Qsel[0][i])
      continue;        /* do not keep    */
    if (op == dp) {
      dp++;
      continue;        /* same particles */
    }
    *dp = *op;         /* Output = Input */
    dp++;
  }
  return 1;
}
/* -------------------------------------------------------------- *\ 
|* copy_onedim_int :
|* Copy selected particles from an one dimensionnal integer array
\* -------------------------------------------------------------- */
int copy_onedim_int(Qsel,oneptr,nbody)
     bool * Qsel[];
     int * oneptr;
     int nbody;
{
  int * op, * dp;
  int  i;

  for (i=0, op=oneptr, dp = oneptr; i<nbody; i++, op++) {
    if (!Qsel[0][i])
      continue;        /* do not keep    */
    if (op == dp) {
      dp++;
      continue;        /* same particles */
    }
    *dp = *op;         /* Output = Input */
    dp++;
  }
  return 1;
}
/* -------------------------------------------------------------- *\ 
|* copy_twodim :
|* Copy selected particles from an two dimensionnal real array
\* -------------------------------------------------------------- */
int copy_twodim(Qsel,twoptr,nbody,step)
     bool * Qsel[];
     real * twoptr;
     int nbody,step; /* phasespace : step = 6 //  acceleration : step = 3 */
{
  real * op, * dp;
  int  i,j;

  for (i=0, op=twoptr, dp = twoptr; i<nbody; i++, op+=step) {
    if (!Qsel[0][i])
      continue;        /* do not keep                       */
    for (j=0;j<step;j++)
      if (op[j]!=dp[j])
	dp[j] = op[j]; /* Output = Input for x,y,z,vx,vy,vz */
    dp+=step;
  }
  return 1;
}

/* -------------------------------------------------------------- *\ 
|* Main program
\* -------------------------------------------------------------- */ 
void nemo_main()         
{
  stream instr,      /* input stream */
    outstr;          /* output stream */

  string timu,                 
    infile,
    headline,
    select_pts[1];

  int coordsys = CSCode(Cartesian, NDIM, 2);

  int    nret[1];

  double * timeptr  = NULL;
  real   * phaseptr = NULL;
  real   * posptr   = NULL;
  real   * velptr   = NULL;
  real   * massptr  = NULL;
  real   * potptr   = NULL;
  real   * accptr   = NULL;
  int    * keyptr   = NULL;

  int    i, j,nbody,nbody_out=0,step,n_step;

  bool   first = TRUE, out=FALSE,
    * Qsel[1];           
    
  int * select_i[1];     
  /*bool within(double, string, double);*/
  /*IECK: within(real,string,real) !!! */

  /* get input parameters */
  infile=getparam("in");
  instr = stropen(getparam("in"), "r");
  select_pts[0] = getparam("select");
   
  timu = getparam("times");
  step = getiparam("step");
  n_step = step;

  outstr = stropen(getparam("out"),"w");
        
  Qsel[0] = NULL;

  get_history(instr);
  put_history(outstr);

  for (;;) {
    /* infinite loop, broken only when ran out of snapshots */
    get_history(instr);                    /* read history */
    while (get_tag_ok(instr,HeadlineTag))
      headline = get_string(instr,HeadlineTag);
       
    if (!get_tag_ok(instr, SnapShotTag)) 
      break; /* check if done */

    get_set(instr, SnapShotTag);

    get_set(instr, ParametersTag);
    if (get_tag_ok(instr,TimeTag)) {
      if (timeptr == NULL)
	timeptr = (double *) allocate(sizeof(double));
      if (timeptr == NULL)
	error("%s: pas assez de memoire\n", getargv0());
      get_data_coerced(instr, TimeTag, DoubleType, timeptr, 0);
      dprintf(1,"Time read : %3.10f\n",*timeptr);
    }

    if (get_tag_ok(instr, NobjTag)) {
      get_data(instr, NobjTag, IntType, &nbody, 0);    
    }
    get_tes(instr, ParametersTag);

    if (get_tag_ok(instr,ParticlesTag) &&
	(timeptr == NULL || streq(timu, "all") ||
	 within((real) *timeptr, timu, (real) TIMEFUZZ))) {

      if (n_step == step) {
	n_step= 1;
                 
	fprintf(stderr,"Copying time steps [%f]...\n",*timeptr);

	/* insertion du select_i pour les champs de saisie */
	for (j=0; j<1; j++)
	  if (Qsel[j]==NULL) {
                       
	    Qsel[j]   = (bool *) allocate (nbody*sizeof(bool));
	    select_i[j] = (int *) allocate (nbody*sizeof(int));
	    if (Qsel[j] == NULL || select_i[j] == NULL)
	      error("%s: not enuf memory\n", getargv0());
                        
	    for (i=0 ; i < nbody ; i++)
	      Qsel[j][i] = FALSE;

	    if (!streq("all",select_pts[j])) {
	      for (i=0; i<nbody; i++) {
		Qsel[j][i] = FALSE;
		select_i[j][i] =-1;
	      }
	      nret[j] = nemoinpi(select_pts[j],
				 select_i[j],nbody);
	      for (i=0; i < nret[j]; i++) {
		Qsel[j][select_i[j][i]]=TRUE;
	      }
	    }
	    else {
	      for (i=0; i<nbody; i++) {
		Qsel[j][i] = TRUE;   
	      }
	      nret[j] = nbody;
	    }
	    free(select_i[j]);
	  } 
	/* Save ParametersTag */
	nbody_out = nret[0];
	out = TRUE;
	save_parameter(outstr,*timeptr,nbody_out);
	get_set(instr, ParticlesTag);
	put_set(outstr, ParticlesTag);

	/* save system coordinates */
	put_data(outstr, CoordSystemTag, IntType, &coordsys, 0);

	/* get masses */
	if (get_tag_ok(instr, MassTag)) {
	  if (massptr == NULL)
	    massptr = (real *) malloc(sizeof(real) * nbody);
                        
	  if (massptr == NULL)
	    error("%s: pas assez de memoire\n", getargv0());
	  get_data_coerced(instr, MassTag, RealType, massptr,
			   nbody, 0);

	  /* keep selectd masses */
	  copy_onedim(Qsel,massptr,nbody);

	  /* save masses */
	  put_data(outstr,MassTag,RealType,massptr,nbody_out,0);
                       
	  /* free allocated memory */
	  free(massptr);
	  massptr = NULL;
	}
	/* get PhaseSpace coordinates */
	if (get_tag_ok(instr, PhaseSpaceTag)) {
	  if (phaseptr == NULL)
	    phaseptr = (real *) allocate(sizeof(real)*2*NDIM*
					 nbody);
                         
	  if (phaseptr == NULL)
	    error("%s: pas assez de memoire\n", getargv0());
	  get_data_coerced(instr, PhaseSpaceTag,
			   RealType,phaseptr,nbody, 2, NDIM, 0);

	  /* keep selected PhaseSpace */
	  copy_twodim(Qsel,phaseptr,nbody,6);

	  /* save phasespace */
	  put_data(outstr, PhaseSpaceTag,
		   RealType,phaseptr,nbody_out, 2, NDIM, 0);

	  /* free allocated memory */
	  free(phaseptr);
	  phaseptr = NULL;
	}
	/* get positions */
	if (get_tag_ok(instr, PosTag)) {
	  if (posptr == NULL)
	    posptr = (real *) malloc(sizeof(real) * nbody * 3);
                          
	  if (posptr == NULL)
	    error("%s: pas assez de memoire\n", getargv0());
	  get_data_coerced(instr, PosTag, RealType, posptr,
			   nbody, 3,0);

	  /* keep selected positions */
	  copy_twodim(Qsel,posptr,nbody,3);

	  /* save positions  */
	  put_data(outstr,PosTag,
		   RealType,posptr,nbody_out, 3, 0);
                       
	  /* free allocated memory */
	  free(posptr);
	  posptr = NULL;
	}
	/* get velocities */
	if (get_tag_ok(instr, VelTag)) {
	  if (velptr == NULL)
	    velptr = (real *) malloc(sizeof(real) * nbody * 3);
                          
	  if (velptr == NULL)
	    error("%s: pas assez de memoire\n", getargv0());
	  get_data_coerced(instr, VelTag, RealType, velptr,
			   nbody, 3,0);

	  /* keep selected velocities */
	  copy_twodim(Qsel,velptr,nbody,3);

	  /* save velocities */
	  put_data(outstr,VelTag,
		   RealType,velptr,nbody_out, 3, 0);
                       
	  /* free allocated memory */
	  free(velptr);
	  velptr = NULL;
	}
	/* get potential */
	if (get_tag_ok(instr, PotentialTag)) {
	  if (potptr == NULL)
	    potptr = (real *) malloc(sizeof(real) * nbody);
                          
	  if (potptr == NULL)
	    error("%s: pas assez de memoire\n", getargv0());
	  get_data_coerced(instr, PotentialTag, RealType, potptr,
			   nbody, 0);

	  /* keep selected potential */
	  copy_onedim(Qsel,potptr,nbody);

	  /* save potential */
	  put_data(outstr,PotentialTag,
		   RealType,potptr,nbody_out,0);
                       
	  /* free allocated memory */
	  free(potptr);
	  potptr = NULL;
	}
	/* get acceleration */
	if (get_tag_ok(instr, AccelerationTag)) {
	  if (accptr == NULL)
	    accptr = (real *) malloc(sizeof(real) * nbody * 3);
                          
	  if (accptr == NULL)
	    error("%s: pas assez de memoire\n", getargv0());
	  get_data_coerced(instr, AccelerationTag, RealType, accptr,
			   nbody, 3,0);

	  /* keep selected accelerations */
	  copy_twodim(Qsel,accptr,nbody,3);

	  /* save accelerations */
	  put_data(outstr,AccelerationTag,
		   RealType,accptr,nbody_out, 3, 0);
                       
	  /* free allocated memory */
	  free(accptr);
	  accptr = NULL;
	}
	/* get Keys */
	if (get_tag_ok(instr, KeyTag)) {
	  if (keyptr == NULL)
	    keyptr = (int *) malloc(sizeof(int) * nbody);
                          
	  if (keyptr == NULL)
	    error("%s: pas assez de memoire\n", getargv0());
	  get_data_coerced(instr, KeyTag, IntType, keyptr,
			   nbody, 0);

	  /* keep selected keys */
	  copy_onedim_int(Qsel,keyptr,nbody);

	  /* save keys */
	  put_data(outstr,KeyTag,
		   IntType,keyptr,nbody_out,0);
                       
	  /* free allocated memory */
	  free(keyptr);
	  keyptr = NULL;
	}
	
	first = FALSE;
	get_tes(instr, ParticlesTag);
	put_tes(outstr, ParticlesTag);
      }
      else {
	n_step++;
      }
    }
    get_tes(instr, SnapShotTag);   
    if (out) {
      put_tes(outstr,SnapShotTag);  
    }
    out = FALSE;  
  }
  strclose(instr);    /* close input file  */
  fclose(outstr);     /* close output file */
  if (first) {
    warning("No snapshots processed");
  } 
  else {
  }
}
/* -------------------------------------------------------------- *\ 
|* End of [snapmask2.c]
\* -------------------------------------------------------------- */ 
