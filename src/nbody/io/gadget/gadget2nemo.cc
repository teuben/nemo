// ============================================================================
// Copyright Jean-Charles LAMBERT - 2006
// e-mail:   Jean-Charles.Lambert@oamp.fr
// address:  Dynamique des galaxies
//           Laboratoire d'Astrophysique de Marseille
//           2, place Le Verrier
//           13248 Marseille Cedex 4, France
//           CNRS U.M.R 6110
// ============================================================================
// gadget2nemo.cc
// gadget2nemo is a program to convert GADGET file to NEMO snapshot.
// ============================================================================
// 08-Nov-2006 : v 2.0 (JCL) - added into NEMO cvs
// 18-Sep-2008 : v 2.1 (WD)  - debugged --- (almost) happy gcc 4.3.1
// ============================================================================
//    NOTE: gcc 4.3.1 warns about conversion from string constant to char* 
//          -- this is s NEMO problem to be addressedd elsewhere.
// ============================================================================

#include <iostream>                                   // C++ I/O
#include <fstream>                                    // C++ file I/O
#include <cstdio>
#include <cstdlib>
#include <snapshot/snapshot.h>
#include <assert.h>

// Endianness I/O class
#include "gadget_data_structure.h"
#include "gadget_endian_tools.h"

// external "C" functions (hey, we are speaking c++ !!!)
extern "C" {
#include <nemo.h>                                     // NEMO basics
  int io_nemo(char *, char *,...);
  int load_snapshot(char *fname, int files);
  int allocate_memory(void);
  int reordering(void);
  int compare (const void * , const void * );
}

using namespace std; // prevent writing statment like 'std::cerr'

#define SPH 0
//------------------------------------------------------------------------------
// NEMO parameters
::string defv[] = {  // use `::'string because of 'using namespace std'
  "in=???\n             GADGET input                                  ",
  "out=???\n            NEMO output                                   ",
  "step=1\n             keep one particle over step                  ",
  "comp=mxv\n           component requested to be saved              ",
  "sort=t\n             sort particles by their Ids                   ",
  "files=1\n            #input gadget files                           ",
  "swap=f\n             #swap data from/to little/big endian          ",
  "VERSION=2.0\n        22-Jun-2005  - JCL                            ",
  NULL
};
::string usage="Convert GADGET snapshot to NEMO snapshot";

int    nbody, * N, step, ntotstep, *iobits=NULL;
float * mass=NULL,  tps, * pos=NULL, * vel=NULL;
bool endian_swap;
::string comp;

//------------------------------------------------------------------------------
// Gadget data structure

t_io_header_1 header1;

int     NumPart, Ngas;

t_particle_data *P;

int *Id;

double  Time, Redshift;

//------------------------------------------------------------------------------
//                            GADGET Functions
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
/* this routine loads particle data from Gadget's default
 * binary file format. (A snapshot may be distributed
 * into multiple files.
 */
int load_snapshot(char *fname, int files)
{
  FILE *fd;
  char   buf[200];
  int    i,k,dummy,ntot_withmasses;
  int    n,pc,pc_new
#if SPH
         ,pc_sph
#endif
         ;
  
#define SKIP iog->ioData((char *) &dummy, sizeof(dummy), 1, GadgetEndianTools::READ);

  cerr << "Files = " << files << "\n";

  for(i=0, pc=1; i<files; i++, pc=pc_new) {

    if(files>1)
      sprintf(buf,"%s.%d",fname,i);
    else
      sprintf(buf,"%s",fname);

    if(!(fd=fopen(buf,"r"))) {
	
      fprintf(stderr,"can't open file `%s`\n",buf);
      std::exit(0);
    }
    GadgetEndianTools * iog = new GadgetEndianTools(fd,endian_swap);

    fprintf(stderr,"reading `%s' ...\n",buf); fflush(stdout);

    SKIP;
    //std::cerr << "header adress =" << &header1 <<"\n";
    iog->ioHeader(&header1,GadgetEndianTools::READ);
    SKIP;

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //do_what_you_want();
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    int localNumPart=0;  // match to #particles in the current local file

    if(files==1) { // Only one input file
      
      for(k=0, NumPart=0, ntot_withmasses=0; k<5; k++) {
	NumPart+= header1.npart[k];
	localNumPart += header1.npart[k];
      }
      Ngas= header1.npart[0];
    }
    else {        // several input files
	
      for(k=0, NumPart=0, ntot_withmasses=0; k<5; k++) {
	NumPart+= header1.npartTotal[k];
	localNumPart += header1.npart[k];
      }
      Ngas= header1.npartTotal[0];
    }

    for(k=0, ntot_withmasses=0; k<5; k++) {
      if(header1.mass[k]==0)
	ntot_withmasses+= header1.npart[k];
    }

    if(i==0) {
      if (ntotstep != 1) {
	ntotstep= NumPart/step;
      } else {
	ntotstep=NumPart;
      }
      char * s[5]={"gaz  ","halo ","disk ","bulge","stars" };
      cerr << "--------------------------------------\n";
      cerr << "Total particles number:               \n";
      cerr << "--------------------------------------\n";
      for (int i=0; i<5; i++ ) {
	cerr << s[i] << " = " << header1.npartTotal[i] << "  step=" <<  header1.npartTotal[i]/step<< "\n";
      }
      cerr << "--------------------------------------\n";
      allocate_memory();
    }

    SKIP;
    for(k=0,pc_new=pc;k<6;k++) {
      for(n=0;n<header1.npart[k];n++) {

	iog->ioData((char *) &P[pc_new].Pos[0], sizeof(float), 3,GadgetEndianTools::READ);
	if ( !(n%step) ) pc_new++;
      }
    }
    //std::cerr << "pc_new = " << pc_new << "\n";
    SKIP;

    SKIP;
    for(k=0,pc_new=pc;k<6;k++) {
      for(n=0;n<header1.npart[k];n++) {

	iog->ioData((char *) &P[pc_new].Vel[0], sizeof(float), 3,GadgetEndianTools::READ);
	  
	if ( !(n%step) ) pc_new++;
      }
    }
    //std::cerr << "pc_new = " << pc_new << "\n";
    SKIP;
    

    SKIP;
    for(k=0,pc_new=pc;k<6;k++) {
      for(n=0;n<header1.npart[k];n++) {
	iog->ioData((char *) &Id[pc_new], sizeof(int), 1,GadgetEndianTools::READ);
	P[pc_new].Id = Id[pc_new];
	if ( !(n%step) ) pc_new++;
      }
    }
    //std::cerr << "pc_new = " << pc_new << "\n";
    SKIP;


    if(ntot_withmasses>0)
      SKIP;
    for(k=0, pc_new=pc; k<6; k++) {
      for(n=0;n<header1.npart[k];n++) {
	P[pc_new].Type=k; // original
		
	if(header1.mass[k]==0)
	  iog->ioData((char *) &P[pc_new].Mass, sizeof(float), 1,GadgetEndianTools::READ);
	else
	  P[pc_new].Mass= header1.mass[k];
	if ( !(n%step) ) pc_new++;
      }
    }
    if(ntot_withmasses>0)
      SKIP;
      
#if SPH
    if(header1.npart[0]>0) {
      SKIP;
      for(n=0, pc_sph=pc; n<header1.npart[0];n++) {
	iog->ioData((char *) &P[pc_sph].U, sizeof(float), 1,GadgetEndianTools::READ);
	if ( !(n%step) ) pc_sph++;
      }
      SKIP;

      SKIP;
      for(n=0, pc_sph=pc; n<header1.npart[0];n++) {

	iog->ioData((char *) &P[pc_sph].Rho, sizeof(float), 1,GadgetEndianTools::READ);
	if ( !(n%step) ) pc_sph++;
      }
      SKIP;

      if(header1.flag_cooling) {
	SKIP;
	for(n=0, pc_sph=pc; n<header1.npart[0];n++) {

	  iog->ioData((char *) &P[pc_sph].Ne, sizeof(float), 1,GadgetEndianTools::READ);
	  if ( !(n%step) ) pc_sph++;
	}
	SKIP;
      }
      else
	for(n=0, pc_sph=pc; n<header1.npart[0];n++) {
	  P[pc_sph].Ne= 1.0;
	  if ( !(n%step) ) pc_sph++;
	}
    }
#endif // if SPH

    fclose(fd);
      
  } // >> for(i=0, pc=1; i<files; i++...


  Time= header1.time;
  Redshift= header1.time;

  return 1;
}



//------------------------------------------------------------------------------
/* this routine allocates the memory for the 
 * particle data.
 */
int allocate_memory()
{
  dprintf(1,"allocating memory...\n");

  if(!(P=(t_particle_data *) malloc(ntotstep*sizeof(t_particle_data))))
    {
      fprintf(stderr,"failed to allocate memory.\n");
      std::exit(0);
    }
  
  P--;   /* start with offset 1 */

  
  if(!(Id=(int *) malloc(ntotstep*sizeof(int))))
    {
      fprintf(stderr,"failed to allocate memory.\n");
      std::exit(0);
    }
  
  Id--;   /* start with offset 1 */

  dprintf(1,"allocating memory...done\n");
  return 1;
}




//------------------------------------------------------------------------------
/* This routine brings the particles back into
 * the order of their ID's.
 * NOTE: The routine only works if the ID's cover
 * the range from 1 to NumPart !
 * In other cases, one has to use more general
 * sorting routines.
 */
int reordering(void)
{
  int i;
  int idsource, idsave, dest;
  t_particle_data psave, psource;


  fprintf(stderr,"reordering....\n");

  for(i=1; i<=ntotstep; i++)
    {
      fprintf(stderr,"i=%08d\tId=%d\n",i,Id[i]);
#if 0
      if ( Id[i] != P[i].Type ) {
	cerr << "Id["<< i << "]="<< Id[i] << "\n";
	cerr  << "P["<< i << "]="<< P[i].Type << "\n";
      }
#endif
      if(Id[i] != i)
	{
	  psource= P[i];
	  idsource=Id[i];
	  dest=Id[i];

	  do
	    {
	      psave= P[dest];
	      idsave=Id[dest];

	      P[dest]= psource;
	      Id[dest]= idsource;
	      fprintf(stderr,"%d\n",i);
	      if(dest == i) 
		break;

	      psource= psave;
	      idsource=idsave;

	      dest=idsource;
	    }
	  while(1);
	  
	}
    }

  fprintf(stderr,"done.\n");

  Id++;   
  free(Id);

  fprintf(stderr,"space for particle ID freed\n");
  return 1;
}

int compare(const void * a, const void * b)
{
  t_particle_data 
    * pa = ( t_particle_data * ) a, 
    * pb = ( t_particle_data * ) b;

  //return (pa->Type - pb->Type);  // sorting according component type
  return (pa->Id - pb->Id);
}
//
//------------------------------------------------------------------------------
// fill_io_nemo
int fill_io_nemo()
{
  static bool first=true;

  if ( first ) {
    first = false;
    // compute nbody
    nbody=ntotstep;
    //for (int i=0; i<6; i++) {
    //  nbody+=header1.npart[i];
    //}
    
    iobits = (int *) allocate(sizeof(int));
    *iobits = TimeBit;
    fprintf(stderr,"NBODY = %d\n",nbody);

    // allocate memory
    if (strchr(comp,'x')) {
      (*iobits) |= PosBit;
      pos  = (float * ) malloc(sizeof(float) * nbody * 3);
      if ( ! pos ) {
	fprintf(stderr,"Unable to allocate memory...\n");
	std::exit(1);
      }
    }
    if (strchr(comp,'v')) {
      (*iobits) |= VelBit;
      vel  = (float * ) malloc(sizeof(float) * nbody * 3);
      if ( ! vel ) {
	fprintf(stderr,"Unable to allocate memory...\n");
	std::exit(1);
      }
    }

    if (strchr(comp,'m')) {
      (*iobits) |= MassBit;
      mass = (float * ) malloc(sizeof(float) * nbody);
      if ( ! mass ) {
	fprintf(stderr,"Unable to allocate memory...\n");
	std::exit(1);
      }
    }
  }
  // read positions/velocities/mass
  for(int k=0,index=0,pc_new=1;k<6;k++) {
    for(int n=0;n<header1.npartTotal[k]/step;n++) {
      assert(index<nbody);

      if ((*iobits) & PosBit) {
	pos[index*3  ]=P[pc_new].Pos[0];
	pos[index*3+1]=P[pc_new].Pos[1];
	pos[index*3+2]=P[pc_new].Pos[2];
      }
      if ((*iobits) & VelBit) {
	vel[index*3  ]=P[pc_new].Vel[0];
	vel[index*3+1]=P[pc_new].Vel[1];
	vel[index*3+2]=P[pc_new].Vel[2];
      }
      if ((*iobits) & MassBit) {
	mass[index   ]=P[pc_new].Mass;
      }

      pc_new++;
      //std::cerr << "My ID = " << P[pc_new].Id << "\n";
      index++;
    }
  }

  return 1;
}
//-----------------------------------------------------------------------------

//------------------------------------------------------------------------------
// main
int main(int argc, char ** argv )
{
  ::string in,out;
  bool sort_part;
  int files;
  //   start  NEMO
  initparam(argv,defv);

  // Get parameters
  in          = getparam("in");
  out         = getparam("out");
  step        = getiparam("step");
  comp        = getparam("comp");
  sort_part   = getbparam("sort");
  files       = getiparam("files");
  endian_swap = getbparam( "swap"  );  
  //cerr << "In = [" << in << "]\n";

  // Read GADGET snapshot
  load_snapshot(in,files);
  //reordering();  /* call this routine only if your ID's are set properly */
  if (sort_part) {
    std::cerr << "Sorting particles by their Ids, may take a while.....\n";
    qsort(P+1,NumPart/step,sizeof(t_particle_data), compare);
  }
  // Fill io_nemo array
  fill_io_nemo();
  tps=0.0;

  // save nemo spnashot
  N = (int * ) malloc(sizeof(int));
  float * tps = (float *) malloc(sizeof(float));
  * tps = header1.time;
  *N = nbody;
  if ( ! io_nemo(out,"save,float,n,t,m,x,v,b,info",
		 &N,&tps,&mass,&pos,&vel,&iobits)) {
    cerr << "Unable to save snapshot [" << out << "]\n";
    cerr << "Aborted\n";
  } else {
  }
  //   finish NEMO
  finiparam();
}
// ----------- End Of [gadget2nemo.cc] ------------------------------------
