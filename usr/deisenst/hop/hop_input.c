/* HOP_INPUT.C, Daniel Eisenstein, 1997 */
/* Based on a paper by Daniel Eisenstein & Piet Hut, 
"HOP: A New Group-Finding Algorithm for N-body Simulations."
See the included documentation or view it at 
http://www.sns.ias.edu/~eisenste/hop/hop_doc.html */

/* Version 1.0 (12/15/97) -- Original Release */
/* Version 1.1 (04/02/01) -- No changes to this file.
			     Tiny bug fix in regroup.c */

/* The routine ReadSimulationFile() is just a wrapper for whatever
routine you need to get your simulation format read from disk and
put in the KD format.  I've included three examples, ReadSimple(),
ReadASCII(), and ReadTPM(), but you can do what you like. */

/* Since you will probably need to write a new version of this,
here's what has to happen:

You are reading from the file fp (which might be stdin and need not be
opened or closed) and putting data into the structure kd.

1) Set kd->nActive:

    kd->nActive = The number of particles you intend to run the algorithm on.

2) Initialize space for kd->p, an array to hold the information for all the
particles.  Do this by:

    kd->p = (PARTICLE *)malloc(kd->nActive*sizeof(PARTICLE));

3) Read in the particle information and put it in the kd->p array.
Note that the array is zero-offset, i.e. it runs from 0 to nActive-1.

You only need to set the position and mass of the particles.  The positions
are at:
	kd->p[j].r[0-2] -- The position information of particle j along
				Cartesian directions 0, 1, and 2.

The masses depend on whether you are compiling with the DIFFERENT_MASSES
flag.  If yes, then kd->p[j].fMass holds the mass of particle j.  If no,
then all the particles are assumed to have the same mass, equal to kd->fMass

The mass and length scales you chose are up to you, but remember that
all your density thresholds will be in those units. 

That's it! */

/* I've included two routines f77read() and f77write() at the bottom
of this file, in case you want to read or write FORTRAN's "unformatted" 
output. */

/* If you only want to run on a subset of the particles, you need to 
adjudicate that within this subroutine and make sure that kd->nActive
and kd->p contain only the particles that you want to include in the 
grouping. */

/* The following variables in the structure KD aren't used by the HOP
program.  You only need to set them if you want to use them for your
own custom purposes (e.g. changing the output to the TIPSY format):
    kd-> nDark, nGas, nStar, nParticles, fTime, bDark, bGas, bStar */

/* Sorry, but I haven't tested ReadASCII or ReadSimple since my files
aren't in that format.  They look ok by eye, but that's been known
to fail :-).  In any case, the point is to give an example of what 
has to be done. */

/* ================================================================ */
/* ====================== ReadSimulationFile() =================== */
/* ================================================================ */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "kd.h"

int ReadTPM(KD kd,FILE *fp);
int ReadASCII(KD kd,FILE *fp);
int ReadSimple(KD kd,FILE *fp);


int ReadSimulationFile(KD kd, FILE *fp)
/* Alter this as needed as described above */
{
    ReadASCII(kd, fp);
    return 0;
}

/* ================================================================ */
/* =================== An Example: ReadSimple() =================== */
/* ================================================================ */

/* Read the following simple format -- unformatted FORTRAN output 
written by the following statements:

	int*4 n_particles
	real*4 pos_x(n_particles), pos_y(n_particles), pos_z(n_particles)
	write(*) n_particles
	write(*) (pos_x(j),j=1,n_particles)
	write(*) (pos_y(j),j=1,n_particles)
	write(*) (pos_z(j),j=1,n_particles)

and that all particles have equal masses, chosen to be 1/n_particles */

int ReadSimple(KD kd,FILE *fp)
{
    int f77read(FILE *f, void *p, int maxbytes);
    int j;
    int header[100]; 	/* Assuming that the first FORTRAN block 
			is smaller than this */
    float *readlist;

    /* First, find out how many particles are involved */
    f77read(fp,header,400);  
    kd->nActive = header[0];  /* The number of particles */

    /* Allocate space to hold their positions */
    kd->p = (PARTICLE *)malloc(kd->nActive*sizeof(PARTICLE));
    assert(kd->p != NULL);

    /* Allocate temporary space to hold the data array */
    readlist = (float *)malloc(4*kd->nActive);

    /* Now read the X positions and transcribe them */
    f77read(fp,readlist,kd->nActive*4);  
    for (j=0; j<kd->nActive; j++)
	kd->p[j].r[0] = readlist[j];	/* Note zero-offset! */

    /* Repeat for Y and Z */
    f77read(fp,readlist,kd->nActive*4);  
    for (j=0; j<kd->nActive; j++)
	kd->p[j].r[1] = readlist[j];	
    f77read(fp,readlist,kd->nActive*4);  
    for (j=0; j<kd->nActive; j++)
	kd->p[j].r[2] = readlist[j];	

    /* Assume the particle mass is 1/kd->nActive */
#ifdef DIFFERENT_MASSES
    for (j=0;j<kd->nActive;j++) kd->p[j].fMass= 1.0/kd->nActive;
#else
    kd->fMass = 1.0/kd->nActive;	
#endif

    /* Give up the temp space */
    free(readlist);
    return kd->nActive;
}

/* ================================================================ */
/* =================== An Example: ReadASCII() =================== */
/* ================================================================ */

/* Read the following format -- an ASCII file with each particle's
information listed line by line:

	Line 1:		N_particles
	Line 2 to N+1:	n X Y Z Mass	

where n is the number of the particle, (X, Y, Z) is the position vector,
and Mass is the mass. */

int ReadASCII(KD kd,FILE *fp)
{
    int j, npart, dummy;
    float pos[3], mass;
    char line[200];	/* We'll read the file line-by-line */
    void f77error(char *s);  /* Report and die */

#ifndef DIFFERENT_MASSES
    /* Our format calls for individual masses, yet we have nowhere to put 
	them and no logical fallback. See the Makefile to compile with
	-DDIFFERENT_MASSES */
    fprintf(stderr,"Don't know what to do with masses.");
    exit(1);
#endif

    /* First, find out how many particles are involved */
    if (fgets(line,200,fp)==NULL) f77error("Unexpected EOF in first line.");
    if (sscanf(line,"%d",&npart)!=1) f77error("Couldn't parse first line.");
    kd->nActive = npart;  /* The number of particles */

    /* Allocate space to hold their positions */
    kd->p = (PARTICLE *)malloc(kd->nActive*sizeof(PARTICLE));
    assert(kd->p != NULL);

    for (j=0;j<npart;j++) {	/* Look at each line */
	if (fgets(line,200,fp)==NULL) f77error("Unexpected EOF.");
	if (sscanf(line,"%d %f %f %f %f", &dummy, pos, pos+1, pos+2, &mass)!=5){
		fprintf(stderr,"Couldn't parse line %d.\n",j+1);
		exit(1);
	}
	/* I won't compare dummy to anything, although this could be a check */
	kd->p[j].r[0] = pos[0];
	kd->p[j].r[1] = pos[1];
	kd->p[j].r[2] = pos[2];
#ifdef DIFFERENT_MASSES 	
	kd->p[j].fMass = mass;
#endif
    }
    return kd->nActive;
}

/* ================================================================ */
/* ====================== An Example: ReadTPM() =================== */
/* ================================================================ */

/* We need to read from Guohong Xu's TPM format */
/* To give info to the user: INFORM("info"); */
#define INFORM(string) printf(string); fflush(stdout)

int ReadTPM(KD kd,FILE *fp)
{
    int f77read(FILE *f, void *p, int maxbytes);
    int header[100];
    int i, j, bl, blocksize;
    float *readlist, masspart;

    f77read(fp,header,400);  /* All the cosmological information */
    f77read(fp,header,8);	/* The particle and block count */
    kd->nActive = header[0];  /* The number of particles */

    /* We're going to use all the particles */
    /* We won't set the following variables; they aren't used unless
    you want to *output* in tipsy format, which isn't my convention: */
    /* kd-> nDark, nGas, nStar, nParticles, fTime, bDark, bGas, bStar */

    kd->p = (PARTICLE *)malloc(kd->nActive*sizeof(PARTICLE));
    assert(kd->p != NULL);

    /* This file format has divided the particle information into 
    header[1] sets.  Each set contains px[], py[], pz[], vx[], vy[],
    and vz[] vectors for the given fraction of the particles. */

    blocksize = header[0]/header[1];
    readlist = (float *)malloc((size_t)(4*blocksize)); 
    assert(readlist!=NULL);
    /* readlist is zero offset for particles */
    
    printf("nActive = %d, blocksize = %d\n", 
	kd->nActive, blocksize);
    INFORM("Reading particles..");
    for (bl=0;bl<header[1];bl++) {
	f77read(fp,readlist,blocksize*4);  /* position_x */
	for (j=0,i=bl*blocksize; j<blocksize; j++,i++)
	    kd->p[i].r[0] = readlist[j];
	f77read(fp,readlist,blocksize*4);  /* position_y */
	for (j=0,i=bl*blocksize; j<blocksize; j++,i++)
	    kd->p[i].r[1] = readlist[j];
	f77read(fp,readlist,blocksize*4);  /* position_z */
	for (j=0,i=bl*blocksize; j<blocksize; j++,i++)
	    kd->p[i].r[2] = readlist[j];
	f77read(fp,readlist,blocksize*4);  /* velocity_x */
	f77read(fp,readlist,blocksize*4);  /* velocity_y */
	f77read(fp,readlist,blocksize*4);  /* velocity_z */
	INFORM(".");
    }
    free(readlist); 

    masspart = 1.0/kd->nActive;	/* All particles have the same mass,
			    chosen so that the average density is 1. */

#ifdef DIFFERENT_MASSES
    for (i=0;i<kd->nActive;i++) kd->p[i].fMass=masspart;
#else
    kd->fMass = masspart;	
#endif

    INFORM("Done!\n");
    return kd->nActive;
}

/* ================================================================ */
/* ===================== Some FORTRAN utilities =================== */
/* ================================================================ */

void f77error(char *s)
{
    fprintf(stderr,"%s\n",s); exit(1);
}

int f77read(FILE *f, void *p, int maxbytes)
/* Read a FORTRAN style block from the given file */
/* maxbytes is the amount of space (in bytes) the pointer p points to */
/* Space must be allocated to read the whole block into p */
/* Return amount read, scream if there's a problem */
/* Reading is done ZERO-OFFSET */
{
    int size, size2;
    if (fread(&size,4,1,f)!=1) 
        f77error("f77read(): Error reading begin delimiter.");
    if (size>maxbytes) 
        f77error("f77read(): Block delimiter exceeds size of storage.");
    if (size<maxbytes) 
        fprintf(stderr,"f77read(): Block size is smaller than size of storage.");
    if (fread(p,1,size,f)!=size) f77error("f77read(): Error reading data.");
    if (fread(&size2,4,1,f)!=1) 
        f77error("f77read(): Error reading end delimiter.");
    if (size!=size2) 
	f77error("f77read(): Delimiters do not match.");
    return size;
}

/*  For completeness.... */
int f77write(FILE *f, void *p, int len)
/* len is number of bytes to be written from p[0..len-1] */
/* Return 0 if successful, 1 if not */
{
    if (fwrite(&len,4,1,f)!=1) return 1;
    if (fwrite(p,1,len,f)!=len) return 1;
    if (fwrite(&len,4,1,f)!=1) return 1;
    return 0;
}

