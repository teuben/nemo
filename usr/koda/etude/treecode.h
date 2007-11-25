/****************************************************************************/
/* TREECODE.H: define various things for treecode.c and treeio.c.           */
/* Copyright (c) 1999 by Joshua E. Barnes, Tokyo, JAPAN.                    */
/*    revised for SPH calculation by Jin Koda, Tokyo, JAPAN. 2000           */
/****************************************************************************/

#ifndef _treecode_h
#define _treecode_h

#include "treedefs.h"

/*
 * Parameters, state variables, and diagnostics for N-body integration.
 */

global string infile;                   /* file name for snapshot input     */
global string outfile;                  /* file name for snapshot output    */
global string savefile;                 /* file name for state output       */
global real freq;                       /* basic integration frequency      */
global real freqout;                    /* data output frequency            */
global real tstop;                      /* time to stop calculation         */
global string headline;                 /* message describing calculation   */
global real tnow;                       /* current value of time            */
global real tout;                       /* time of next output              */
global int nstep;                       /* number of time-steps             */
global int nbody;                       /* number of bodies in system       */
global bodyptr bodytab;                 /* points to array of bodies        */

global real rscale;                     /* scale of radius [kpc]            */
global real mscale;                     /* scale of mass [Msun]             */
global real tscale;                     /* scale of time [yr]               */
global real massdk;                     /* total disk mass in system unit   */
global real dtcfl;                      /* t.step crit. for CFL             */
global real dtvel;                      /* t.step crit. for velocity        */
global real dtacc;                      /* t.step crit. for acceleration    */

/*
 * Prototypes for I/O routines.
 */

void inputdata(void);                   /* read initial data file           */
void startoutput(void);                 /* open files for output            */
void forcereport(void);                 /* report on force calculation      */
void output(void);                      /* perform output operation         */
void savestate(string);                 /* save system state                */
void restorestate(string);              /* restore system state             */

#endif /* ! _treecode_h */
