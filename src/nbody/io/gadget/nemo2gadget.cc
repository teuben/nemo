// -----------------------------------------------------------------------------
// nemo2gadget.cc
// -----------------------------------------------------------------------------
#include <iostream>                                   // C++ I/O
#include <fstream>                                    // C++ file I/O
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>

#include "gadget_data_structure.h"
#include "gadget_endian_tools.h"

extern "C" {
#include <nemo.h>                                     // NEMO basics
  int io_nemo(char *, char *,...);
  void savepositions_ioformat1(char * buf);
  size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE *stream);
  int get_chunk(int, int ,int );
}

using namespace std; // prevent writing statment like 'std::cerr'

#define CHUNK(comp,maxfile,current) (((current)+1 == (maxfile)) ? ( (comp)-((comp)/current)): ((comp)/current))
//------------------------------------------------------------------------------
// NEMO parameters
::string defv[] = {  // use `::'string because of 'using namespace std'
  "in=???\n          NEMO   input                                 ",
  "out=???\n         GADGET output                                ",
  "ng=0\n            #gas     particles                           ",
  "nh=0\n            #halo    particles                           ",
  "nd=0\n            #disk    particles                           ",
  "nb=0\n            #bulge   particles                           ",
  "ns=0\n            #stars   particles                           ",
  "files=1\n         #output gadget files                         ",
  "swap=f\n          #swap data from/to little/big endian         ",
  "VERSION=1.1\n     22-Jan-2004  - JCL                           ",
  NULL
};
::string usage="Convert NEMO snapshot to GADGET snapshot";

int    * nbody=NULL;
float * mass=NULL,  * tps = NULL, * pos=NULL, * vel=NULL;
int ng,nh,nd,nb,ns;
int files;
bool endian_swap;

//------------------------------------------------------------------------------
// Gadget data structure

t_io_header_1 header1;

int     NumPart, Ngas;

t_particle_data *P;

int *Id;

double  Time, Redshift;

//#define MAX(A,B) ((A)>(B)?(A):(B))
int get_chunk(int comp, int maxfile, int current) 
{
  int ret;
  int dd= comp/maxfile;

  if (current+1 == maxfile) {
    if ((comp-dd*maxfile) != 0 ) {
      ret = dd+comp-dd*maxfile;
    }
    else {
      ret = dd;
    }
  }
  else {
    ret = dd;
  }
  return ret;
}
//------------------------------------------------------------------------------
//                            GADGET Functions
//------------------------------------------------------------------------------

/* This function writes a snapshot of the particle ditribution to
 * one file using Gadget's default file format.
 * Each snapshot file contains a header first, then particle positions, 
 * velocities and ID's.
 * Then particle masses are written for those particle types with zero entry in
 * MassTable.
 * After that, first the internal energies u, and then the density is written
 * for the SPH particles.
 * Finally, if cooling is enabled, the mean molecular weight is written for the gas
 * particles. 
 */
void savepositions_ioformat1(char * buf)
{
  FILE *fd;
  float dummy[3];
  int i,k;
  int   blklen,masscount;
  double a3inv;
#ifdef COOLING
  double ne, nh0;
#endif

#define ORI_IO 0

#if ORI_IO
# define BLKLEN my_fwrite(&blklen, sizeof(blklen), 1, fd);
#else
# define BLKLEN iog->ioData((char *) &blklen, sizeof(blklen), 1, GadgetEndianTools::WRITE);
#endif

  //  if(All.ComovingIntegrationOn)
  //  a3inv=  1/(All.Time*All.Time*All.Time);
  //else
  a3inv=  1.0;

  //NumPart = *nbody;
  
  // MaxParticles
  header1.npartTotal[0] = ng;
  header1.npartTotal[1] = nh;
  header1.npartTotal[2] = nd;
  header1.npartTotal[3] = nb;
  header1.npartTotal[4] = ns;
  header1.npartTotal[5] = 0;

  int offset[6];

  offset[0] = 0;
  offset[1] = ng;
  offset[2] = ng + nh;
  offset[3] = ng + nh + nd;
  offset[4] = ng + nh + nd + ns;
  offset[5] = 0;

  for (int nf=0; nf<files; nf++) {
    char tmp[200];
    GadgetEndianTools * iog;
    if (files==1) {
      sprintf(tmp,"%s",buf);
      cerr << "Processing file [" << tmp << "]\n";
    } 
    else {
      sprintf(tmp,"%s.%d",buf,nf);
      cerr << "Processing file [" << tmp << "]\n";
    }
    if((fd=fopen(tmp,"w")))
      {
	iog = new GadgetEndianTools(fd,endian_swap);
	header1.npart[0]= get_chunk(ng,files,nf);
	//cerr << ">>" << header1.npart[0] <<" \n";
	header1.npart[1]= get_chunk(nh,files,nf);
	//cerr << ">>" << header1.npart[1] <<" \n";
	header1.npart[2]= get_chunk(nd,files,nf);
	//cerr << ">>" << header1.npart[2] <<" \n";
	header1.npart[3]= get_chunk(nb,files,nf);
	//cerr << ">>" << header1.npart[3] <<" \n";
	header1.npart[4]= get_chunk(ns,files,nf);
	//cerr << ">>" << header1.npart[4] <<" \n";
	header1.npart[5]= 0;

	NumPart = 0;
	for(i=0;i<6;i++) {
	  header1.mass[i]=0;
	  NumPart+=header1.npart[i];
	}

	for(i=0, masscount=0; i<5; i++)
	  {
	    masscount+= header1.npart[i];
	  }

	header1.time= *tps;

	//  if(All.ComovingIntegrationOn)
	//	header1.redshift=1.0/All.Time - 1.0;
	//else
	header1.redshift=0;  

      
	header1.flag_sfr=0;
	header1.flag_feedback=0;
	header1.flag_cooling= 0;
#ifdef COOLING
	header1.flag_cooling= 1;
#endif
	header1.num_files= files;
	header1.BoxSize= 0.;
	header1.Omega0=  0.;
	header1.OmegaLambda= 0.;
	header1.HubbleParam= 0.;
      
	blklen=sizeof(header1);
	BLKLEN;
#if ORI_IO
	my_fwrite(&header1, sizeof(header1), 1, fd);
#else
	iog->ioHeader(&header1,GadgetEndianTools::WRITE);
#endif
	BLKLEN;


	blklen=NumPart*3*sizeof(float);

	BLKLEN;
	//for(i=0;i<NumPart;i++)
	for (int j=0; j< 6; j++) {                  // loop on all components
	  for (int i=0; i<header1.npart[j]; i++) {  // loop on component's particles
	    for(k=0;k<3;k++)                        // 3dim
	      dummy[k]=pos[(offset[j]+i)*3+k];
#if ORI_IO
	    my_fwrite(dummy,sizeof(float),3,fd);
#else

	    iog->ioData((char *) dummy,sizeof(float),3,GadgetEndianTools::WRITE);
#endif
	  }
	}
	BLKLEN;


	BLKLEN;
	//for(i=0;i<NumPart;i++)
	for (int j=0; j< 6; j++) {                  // loop on all components
	  for (int i=0; i<header1.npart[j]; i++) {  // loop on component's particles
	    for(k=0;k<3;k++)                        // 3dim
	      dummy[k]=vel[(offset[j]+i)*3+k];
#if ORI_IO
	    my_fwrite(dummy,sizeof(float),3,fd);
#else
	    iog->ioData((char *) dummy,sizeof(float),3,GadgetEndianTools::WRITE);
#endif
	  }
	}
	BLKLEN;
 

	blklen=NumPart*sizeof(int);
	BLKLEN;
	int myid=1;
	for (int k=0; k< 6; k++) {
	  for(int n=0;n<header1.npart[k];n++) {
#if ORI_IO
	    my_fwrite(&k,sizeof(int),1,fd);
#else
	    iog->ioData((char *) &myid,sizeof(int),1,GadgetEndianTools::WRITE);
	    myid++;
#endif
	  }
	}
	BLKLEN;

	blklen=masscount*sizeof(float);
	if(masscount)
	  BLKLEN;
	//for(i=0;i<NumPart;i++)
	for (int j=0; j< 6; j++) {                  // loop on all components
	  for (int i=0; i<header1.npart[j]; i++) {  // loop on component's particles
	    dummy[0]= mass[(offset[j]+i)];
#if ORI_IO
	    my_fwrite(dummy,sizeof(float),1,fd);
#else
	    iog->ioData((char *) dummy,sizeof(float),1,GadgetEndianTools::WRITE);
#endif
	  }
	}
	if(masscount)
	  BLKLEN;

	if(header1.npart[0])  // gas particles
	  {
	    fprintf(stderr,"Unable to proceed GAS particles : NOT implemented....\n");
	    exit(1);
#if 0
	    blklen=N_gas*sizeof(float);
	    BLKLEN;
	    for(i=1;i<=N_gas;i++)
	      {
		dummy[0]=SphP[i].EgySpecPred;
#if ORI_IO
		my_fwrite(dummy,sizeof(float),1,fd);
#else
		iog->ioData((char *) dummy,sizeof(float),1,GadgetEndianTools::WRITE);
#endif
	      }
	    BLKLEN;


	    blklen=N_gas*sizeof(float);  /* added density  */
	    BLKLEN;
	    for(i=1;i<=N_gas;i++)
	      {
		dummy[0]=SphP[i].DensityPred;
#if ORI_IO
		my_fwrite(dummy,sizeof(float),1,fd);
#else
		iog->ioData((char *) dummy,sizeof(float),1,GadgetEndianTools::WRITE);
#endif
	      }
	    BLKLEN;

#ifdef COOLING
	    blklen=N_gas*sizeof(float);  /* electron abundance */
	    BLKLEN;
	    for(i=1;i<=N_gas;i++)
	      {
		dummy[0]= SphP[i].Ne;
#if ORI_IO
		my_fwrite(dummy,sizeof(float),1,fd);
#else
		iog->ioData((char *) dummy,sizeof(float),1,GadgetEndianTools::WRITE);
#endif
	      }
	    BLKLEN;


	    blklen=N_gas*sizeof(float);  /* neutral hydrogen */
	    BLKLEN;
	    for(i=1;i<=N_gas;i++)
	      {
		ne= SphP[i].Ne;

		AbundanceRatios(SphP[i].EgySpecPred, SphP[i].DensityPred*a3inv,
				&ne, &nh0);
		dummy[0]= nh0;
#if ORI_IO
		my_fwrite(dummy,sizeof(float),1,fd);
#else
		iog->ioData((char *) dummy,sizeof(float),1,GadgetEndianTools::WRITE);
#endif
	      }
	    BLKLEN;
#endif
	    blklen=N_gas*sizeof(float);  /* hsml  */
	    BLKLEN;
	    for(i=1;i<=N_gas;i++)
	      {
		dummy[0]=SphP[i].Hsml;
#if ORI_IO
		my_fwrite(dummy,sizeof(float),1,fd);
#else
		iog->ioData((char *) dummy,sizeof(float),1,GadgetEndianTools::WRITE);
#endif		
	      }
	    BLKLEN;
#endif
	  }
	fclose(fd);
      }
    else
      {
	fprintf(stdout,"Error. Can't write in file '%s'\n", buf);
	//endrun(10);
	exit(1);
      }
    // Offset update
    for (int j=0; j< 6; j++) {                  // loop on all components
      offset[j]+=header1.npart[j];
    }

  } // << for (int nf...)
}


// This catches I/O errors occuring for my_my_fwrite(). In this case we better stop.
size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE *stream)
{
  size_t nwritten;

  if((nwritten=fwrite(ptr, size, nmemb, stream))!=nmemb)
    {
      printf("I/O error (fwrite) on has occured.\n");
      fflush(stdout);
      //endrun(777);
      exit(1);
    }
  return nwritten;
}


//-----------------------------------------------------------------------------

//------------------------------------------------------------------------------
// main
int main(int argc, char ** argv )
{
  ::string in,out;

  //   start  NEMO
  initparam(argv,defv);

  // Get parameters
  in   = getparam ( "in"    );
  out  = getparam ( "out"   );
  ng   = getiparam( "ng"    );
  nh   = getiparam( "nh"    );
  nd   = getiparam( "nd"    );
  nb   = getiparam( "nb"    );
  ns   = getiparam( "ns"    );
  files= getiparam( "files" );
  endian_swap = getbparam( "swap"  );
  // save nemo spnashot
  if ( ! io_nemo(in,"read,float,n,t,m,x,v,info",
		 &nbody,&tps,&mass,&pos,&vel)) {
    cerr << "Unable to read snapshot [" << in << "]\n";
    cerr << "Aborted\n";
  } else {
    savepositions_ioformat1(out);
  }
  //   finish NEMO
  finiparam();
}
// ----------- End Of [nemo2gadget.cc] ------------------------------------
