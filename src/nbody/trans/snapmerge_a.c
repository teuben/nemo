/* -------------------------------------------------------------- *\     
|* snapmerge_a.c :	JC LAMBERT	        	V1.21            
|*                                                                       
|* Summary : Snapmerge_a merge the second input snapshot at the end      
|* of the first snapshot but without making a copy of the first snapshot.
|* This is very convenient when the first snapshot is very large and     
|* prevent problem of space disk.                                        
|*                                                                       
|* the program first go to the end of the first snapshot and get nbody   
|* and time. Then, it reads each timestep of the second snapshot and     
|* starts to append to the first snapshot only if the time of the current
|* timestep is greater than the time of the first, and if the nbody values
|*  are identical.                                                       
|*                                                                       
|* the program print on the stdout the following codes :                 
|* -1 : merging did not occur bc of error.                               
|*  0 : the second snapshot does not have a time step greater than the last
|*      time step of the first snapshot, no merging.                     
|*  1 : merging have been successfully completed.                        
|*                                                                       
|* V1.0  : 10-Jun-96 Created                                             
|* V1.01 : 08-Jul-02 PosTag and VelTag added                             
|* v1.2  : 23-Jan-06 EpsTag added                                        
|* V1.21 : 28-Nov-06 in NEMO cvs                                         
\* -------------------------------------------------------------- */

/* -------------------------------------------------------------- *\
|* Include files                                                   
\* -------------------------------------------------------------- */
#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>		
#include <filestruct.h>
#include <math.h>

#include <snapshot/snapshot.h>	
#include <snapshot/body.h>

/* -------------------------------------------------------------- *\
|* Shared variables                                                 
\* -------------------------------------------------------------- */ 
string defv[] = {
  "in1=???\n			Input NEMO snapshot",
  "in2=???\n			Input NEMO snapshot to append after in1",
  "VERSION=1.22\n               6-Jun-08 pjt",
  NULL,
};

string usage="Merge the second input snapshot at the end of the first one";

#define TIMEFUZZ 0.00001

real      
  * timeptr  = NULL,
  * posptr   = NULL,
  * velptr   = NULL,
  * phaseptr = NULL,
  * massptr  = NULL,
  * potptr   = NULL,
  * accptr   = NULL,
  * epsptr   = NULL;

int      * nbodyptr = NULL;
int        nbody;

/* -------------------------------------------------------------- *\
|* save_parameter :                                                 
\* -------------------------------------------------------------- */
int save_parameter(outstr,timu,nbody)
     stream outstr;
     real * timu;
     int  nbody;
{
  put_set(outstr, SnapShotTag);		/*   start snapshot output  */
  put_set(outstr, ParametersTag);
  put_data(outstr, TimeTag, RealType, timu, 0);
  put_data(outstr, NobjTag, IntType, &nbody, 0);
  put_tes(outstr, ParametersTag);  
}
/* -------------------------------------------------------------- *\
|* append_snap :                                                    
|* Append the snapshot 2 after the snapshot 1.                      
\* -------------------------------------------------------------- */
int append_snap(stream outstr, 
                stream instr,
                int nbody,
                real timu,
                string infile1,
                string infile2)
{ 
  int nbody2;
  bool out_f=FALSE,one=FALSE;
  string headline;

  int coordsys = CSCode(Cartesian, NDIM, 2);
 
  for (;;) {
    /* infinite loop, broken only when ran out of snapshots */
    get_history(instr);                    /* read history */
    while (get_tag_ok(instr,HeadlineTag))
      headline = get_string(instr,HeadlineTag); 
       
    if (!get_tag_ok(instr, SnapShotTag)) 
      break; /* check if done */

    get_set(instr, SnapShotTag);

    get_set(instr, ParametersTag);

    get_data_time(instr, RealType, sizeof(real),(real*) &timeptr); 
    get_data_nbody(instr, IntType, sizeof(int),(int*) &nbodyptr);
    nbody2 = *nbodyptr; 

    if (nbody2 != nbody){
      dprintf(0,"Inconsistent nbody : %s (%d) <> %s (%d)\n",
	      infile1,nbody,infile2,nbody2);
      return -1;

    }
    get_tes(instr, ParametersTag);

    if (*timeptr > timu) {
      out_f=TRUE; one=TRUE;
      save_parameter(outstr,timeptr,nbody);
      
      get_set(instr, ParticlesTag);
      put_set(outstr, ParticlesTag);
      
      /* Save system coordinates */
      put_data(outstr, CoordSystemTag, IntType, &coordsys, 0);
      
      /* recuperation de la masse */
      
      if (!get_data_mass(instr,RealType,*nbodyptr,sizeof(real),
			 &massptr)) {
	dprintf(1,"No Mass \n");
        
      } 
      else {
	/* Save the mass */
	put_data(outstr,MassTag,RealType,massptr,nbody,0);
      }             
      
      /* get phasespace coordinates */
      if (!get_data_phase(instr,RealType,*nbodyptr,sizeof(real), &phaseptr,NDIM)) {
	dprintf(1,"No PhaseSpace !!\n");
	
      }
      else {
	/* Save phasespace coordinates */
	
	put_data(outstr, PhaseSpaceTag,
		 RealType,phaseptr,nbody, 2, NDIM, 0);
      }  
      
      /* getting positions */
      if (!get_data_pos(instr,RealType,*nbodyptr,sizeof(real),
			&posptr,NDIM)) {
	dprintf(1,"No Positions !!\n");
                   
      }
      else {
	/* Save positions */
	put_data(outstr, PosTag,
		 RealType,posptr,nbody, 3, 0);
      }  
      /* getting velocities */
      if (!get_data_vel(instr,RealType,*nbodyptr,sizeof(real),
			&velptr,NDIM)) {
	dprintf(1,"No Velocities !!\n");
                   
      }
      else {
	/* Save velocities */
	put_data(outstr, VelTag,
		 RealType,velptr,nbody, 3, 0);
      }  

      /* recuperation du potentiel */
      if (!get_data_pot(instr,RealType,*nbodyptr,sizeof(real),
			&potptr)) {
      }
      else { /* Save potential */
	put_data(outstr,PotentialTag,
		 RealType,potptr,nbody,0);
      }

      /* recuperation de l'acceleration */
      if (!get_data_acc(instr,RealType,*nbodyptr,sizeof(real),
			&accptr,3)) {
      }
      else { /* Save acceleration */
	put_data(outstr,AccelerationTag,
		 RealType,accptr,nbody,3,0);
      }

      /* try to get EPS fields */
      if (!get_data_eps(instr,RealType,*nbodyptr,sizeof(real),
			&epsptr)) {
      }
      else { /* Save Epsilons */
	put_data(outstr,EpsTag,
		 RealType,epsptr,nbody,0);
      }      

      get_tes(instr, ParticlesTag);
      put_tes(outstr,ParticlesTag);
    }
    else {
      dprintf(0,"Skipping time !!  %s (%f) <= %s (%f)\n",
	      infile2,*timeptr,infile1,timu);
    }
    get_tes(instr, SnapShotTag);
    if (out_f)
      put_tes(outstr, SnapShotTag);
    out_f = FALSE;
  }
  strclose(instr);
  strclose(outstr);

  if (one) 
    return 1;
  else     
    return 0;
    
}
/* -------------------------------------------------------------- *\
|* goto_end_1st_snap :                                              
|* Go to the end of the 1st snapshot, return the last time and the  
* nbody.                                                            
\* -------------------------------------------------------------- */ 
int goto_end_1st_snap(stream   instr,
                      int    * nbody, 
                      real   * timu)
{
  bool first=TRUE;
  string headline;

  for (;;) {
    /* infinite loop, broken only when ran out of snapshots */
    get_history(instr);                    /* read history */
    while (get_tag_ok(instr,HeadlineTag))
      headline = get_string(instr,HeadlineTag); 
       
    if (!get_tag_ok(instr, SnapShotTag)) 
      break; /* check if done */

    get_set(instr, SnapShotTag);

    get_set(instr, ParametersTag);
    get_data_time(instr, RealType, sizeof(real),(real*) &timeptr); 
    get_data_nbody(instr, IntType, sizeof(int),(int*) &nbodyptr);
    *nbody = *nbodyptr;
    *timu  = *timeptr;
    first=FALSE;
    get_tes(instr, ParametersTag);

    get_tes(instr, SnapShotTag);
  }
  strclose(instr);

  dprintf(1,"Nbody = %d, Time = %f\n",*nbody,*timu);

  if (!first)
    return 1;
  else
    return 0;

}
/* -------------------------------------------------------------- *\
|* Main program                                                     
\* -------------------------------------------------------------- */ 
nemo_main()         
{
  stream in1,out,  
    in2;     
  string in_out_file, infile2;
  int    nboby;
  real   timu;
  bool   first = TRUE, compute;

  /* get input parameters */
  in_out_file = getparam("in1");
  infile2     = getparam("in2");

  /* open the input file 1 in reading */
  in1 = stropen(in_out_file,"r");

  /* open the input file 1 in append */
  out = stropen(in_out_file,"a");

  /* open the input file 2 in reading */
  in2 = stropen(infile2,"r");

  /* go to the end of the first snapshot */
  if (!goto_end_1st_snap(in1,&nbody,&timu)) {
    dprintf(0,"No parameters field in <%s>\n",in_out_file);
    printf("-1");
    exit(1);
  }

  /* append the second snapshot after the first */
  printf("%d\n",append_snap(out,in2,nbody,timu,in_out_file,infile2));
} 
/* -------------------------------------------------------------- *\
|* End of snapmerge_a.c                                             
\* -------------------------------------------------------------- */ 
