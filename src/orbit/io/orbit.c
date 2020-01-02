/* orbit.c - write_orbit, read_orbit, allocate_orbit, copy_orbit, list_orbit */

/*------------------------------------------------------------------------------
 * ORBIT      proposed format for orbit paths
 *              documentation in orbit.3
 *
 *  13-Jul-87   V1.0 created                            P.J. Teuben
 *  28-jul-87   V2.0 new orbit(5)                               PJT
 *   3-May-88   ---- added few dprintf()'s and data history     PJT
 *   2-jun-88   ---- newly compiled for new filestruct          PJT
 *                      get/put_hist now
 *  22-may-90   new copy_orbit, debug read_orbit                PJT
 *  11-jun-90   list_orbit output
 *  19-dec-90   read_orbit uses allocate_orbit for allocate()   PJT
 *		  but now there is no method to save allocated space PJT
 *  24-may-92   V3.0 added potential                            PJT
 *  23-nov-93   V3.1 but fixed rather emberassing bug           pjt
 *   6-dec-93	V3.2 added maxsteps				pjt
 *  22-feb-95   V3.2a ansi					pjt
 *  17-apr-95   V3.2b compacted needed header files		pjt
 *  19-apr-96       c initialize strings to "" instead of NULL  pjt
 *   1-mar-03   V3.3  added stuff for iom_err			pjt
 *  25-jul-13   V4.0  added support for key                     pjt
 *  12-nov-2015 V5.0  some support for split Pos/Vel            pjt
 *  10-dec-2019 V5.1  optional support for PHI/ACC              pjt
 *------------------------------------------------------------------------------
 */

#include <stdinc.h>
#include <orbit.h>

#define ISSTRING(s) ( (s)!=NULL && *(s)!=0 )

/*------------------------------------------------------------------------------
 *  WRITE_ORBIT: writes out one orbit
 *------------------------------------------------------------------------------
 */

void write_orbit (stream outstr, orbitptr optr)
{
    /* No history written here, user should do this */
    put_set (outstr,OrbitTag);
        put_set (outstr,ParametersTag);
            put_data (outstr,NdimTag,    IntType, &(Ndim(optr)),   0);
            put_data (outstr,MassTag,    RealType,&(Masso(optr)),   0);
            put_data (outstr,KeyTag,     IntType, &(Key(optr)),   0);
            put_data (outstr,IOMTag,     RealType,  IOM(optr), Ndim(optr), 0);
            put_data (outstr,IOMERRTag,  RealType,  IOMERR(optr), Ndim(optr), 0);
	    /* note we don't write MAXsteps(optr) to disk */
            put_data (outstr,NstepsTag,  IntType, &(Nsteps(optr)), 0);
        put_tes (outstr,ParametersTag);

        put_set (outstr,PotentialTag);
            if (ISSTRING(PotName(optr))) put_string (outstr, PotNameTag, PotName(optr));
            if (ISSTRING(PotPars(optr))) put_string (outstr, PotParsTag, PotPars(optr));
            if (ISSTRING(PotFile(optr))) put_string (outstr, PotFileTag, PotFile(optr));
        put_tes (outstr,PotentialTag);

        put_set (outstr,PathTag);
            put_data (outstr,TimePathTag, RealType, 
                                TimePath(optr), Nsteps(optr), 0);
            put_data (outstr,PhasePathTag,RealType, 
                                PhasePath(optr), Nsteps(optr), 2, Ndim(optr), 0);
#ifdef ORBIT_PHI
	    put_data (outstr,PhiPathTag,RealType,
		      PhiPath(optr), Nsteps(optr), 0);
	    put_data (outstr,AccPathTag,RealType,
		      AccPath(optr), Nsteps(optr), Ndim(optr), 0);
	    
#endif	    
        put_tes (outstr,PathTag);
    put_tes (outstr,OrbitTag);
}

/*------------------------------------------------------------------------------
 * READ_ORBIT: try and read another orbit
 *            returns 0 if not, or an error occurred,  1 if OK
 *------------------------------------------------------------------------------
 */
 
int read_orbit (stream instr, orbitptr *optr)
{
    int  maxsize, ndim, nsteps;
    orbitptr otmp;

    get_history(instr);     /* always scan for history; allow sandwiched */

    // @todo if Orbit not found, continue until one if found or EOF
    if (!get_tag_ok (instr,OrbitTag))
        return 0;                      /* not another orbit available */
        
    if (*optr!=NULL)
        maxsize = Size(*optr);
    else
        maxsize = 0;
        
    get_set (instr,OrbitTag);
        get_set (instr,ParametersTag);
            get_data (instr,NdimTag,    IntType, &ndim,   0);
            get_data (instr,NstepsTag,  IntType, &nsteps, 0);
            if (ndim*nsteps > maxsize) {	      /* allocate more space */
                if (*optr) free_orbit(*optr);         /* free the old one    */
                allocate_orbit(&otmp,ndim,nsteps);
                *optr = otmp;
		Size(*optr) = ndim*nsteps;
		dprintf(1,"Reallocated orbit with %d space\n",ndim*nsteps);
            } else {
                Ndim(*optr) = ndim;
                Nsteps(*optr) = nsteps;
            }
            get_data (instr, MassTag, RealType, &(Masso(*optr)), 0);
	    if (get_tag_ok(instr,KeyTag))
	      get_data (instr, KeyTag, IntType, &Key(*optr), 0);
	    else
	      Key(*optr) = 0;
            get_data (instr, IOMTag, RealType, IOM(*optr), Ndim(*optr), 0);
            if (get_tag_ok(instr,IOMERRTag))
                get_data (instr, IOMERRTag, RealType, IOMERR(*optr), Ndim(*optr), 0);
            else
                IE1(*optr) = IE3(*optr) = IE3(*optr) = 0.0;
        get_tes (instr,ParametersTag);

#if 0
        PotName(*optr) = NULL;
        PotPars(*optr) = NULL;
        PotFile(*optr) = NULL;
#else
        PotName(*optr) = "";
        PotPars(*optr) = "";
        PotFile(*optr) = "";
#endif        
	if (get_tag_ok(instr,PotentialTag)) {
            get_set(instr,PotentialTag);
            if(get_tag_ok(instr,PotNameTag)) 
                PotName(*optr) = get_string(instr,PotNameTag);
            if(get_tag_ok(instr,PotParsTag)) 
                PotPars(*optr) = get_string(instr,PotParsTag);
            if(get_tag_ok(instr,PotFileTag)) 
                PotFile(*optr) = get_string(instr,PotFileTag);
            get_tes(instr,PotentialTag);
	} 

        get_set (instr,PathTag);
            get_data (instr,TimePathTag,    RealType, 
                      TimePath(*optr), Nsteps(*optr), 0);
            get_data (instr,PhasePathTag,RealType, 
                      PhasePath(*optr), Nsteps(*optr), 2, Ndim(*optr), 0);
#ifdef ORBIT_PHI
	    dprintf(1,"Reading orbit Phi/Acc\n");
	    if (get_tag_ok(instr, PhiPathTag)) {
	      get_data (instr,PhiPathTag,    RealType, 
                      PhiPath(*optr), Nsteps(*optr), 0);
	      get_data (instr,AccPathTag,RealType, 
		      AccPath(*optr), Nsteps(*optr), Ndim(*optr), 0);
	    }
#endif	    
        get_tes (instr,PathTag);
    get_tes (instr,OrbitTag);
    return 1;
}

void free_orbit(orbitptr optr)
{
    free( (char *) TimePath(optr) );
    free( (char *) PhasePath(optr) );
    free( (char *) IOM(optr) );
    free( (char *) optr );
}

/*------------------------------------------------------------------------------
 * ALLOCATE_ORBIT: allocate space for an orbit path
 *      returns 0 if failure due to lack of memory or so
 *      returns 1 if all seems OK.
 *------------------------------------------------------------------------------
 */
 
int allocate_orbit(orbitptr *optr, int ndim, int nsteps)
{
    *optr = (orbitptr ) allocate(sizeof(orbit));
    dprintf (1,"Allocate_Orbit @ %x with ndim=%d nsteps=%d\n",
                *optr, ndim, nsteps);
    Ndim(*optr) = ndim;
    IOM(*optr)    = (real *) allocate(ndim*sizeof(real));
    dprintf (2,"   IOM allocated @ %x\n",IOM(*optr)); 
    IOMERR(*optr) = (real *) allocate(ndim*sizeof(real));
    dprintf (2,"   IOMERR allocated @ %x\n",IOMERR(*optr)); 
    Nsteps(*optr) = nsteps;
    TimePath(*optr) = (real *) allocate(nsteps*sizeof(real));
    dprintf (2,"   Timepath allocated @ %x\n",TimePath(*optr));
    PhasePath(*optr) = (real* ) allocate(nsteps*ndim*2*sizeof(real));
    dprintf (2,"  Phasepath allocated @ %x\n",PhasePath(*optr));
#ifdef ORBIT_PHI
    PhiPath(*optr) = (real *) allocate(nsteps*sizeof(real));
    AccPath(*optr) = (real* ) allocate(nsteps*ndim*sizeof(real));
#endif    
    Size(*optr) = ndim*nsteps;
    MAXsteps(*optr) = nsteps;  /* ??? */

    return 1;      /* succes */
}

/*------------------------------------------------------------------------------
 * COPY_ORBIT: copy an orbit, assumes destination has been allocated.
 *------------------------------------------------------------------------------
 */
 
void copy_orbit(orbitptr iptr, orbitptr optr)
{
        int i, nsteps;
        
        if (optr==NULL) error("copy_orbit: Output orbit not allocated");
	
	if (Nsteps(optr) < Nsteps(iptr)) 
	  error("copy_orbit: not enough space in orbit to copy to: out->%d < in->%d",
		Nsteps(optr),Nsteps(iptr));
        
        Ndim(optr) = Ndim(iptr);
        CoordSys(optr) = CoordSys(iptr);
        Masso(optr) = Masso(iptr);
	Key(optr) = Key(iptr);
        IOM(optr) = IOM(iptr);                 /* needs a copy, not a link ?? */
        IOMERR(optr) = IOMERR(iptr);
        nsteps = Nsteps(optr) = Nsteps(iptr);
        for (i=0; i<nsteps; i++) {
            Torb(optr,i) = Torb(iptr,i);
            Xorb(optr,i) = Xorb(iptr,i);
            Yorb(optr,i) = Yorb(iptr,i);
            Zorb(optr,i) = Zorb(iptr,i);
            Uorb(optr,i) = Uorb(iptr,i);
            Vorb(optr,i) = Vorb(iptr,i);
            Worb(optr,i) = Worb(iptr,i);
#ifdef ORBIT_PHI	    
	    Porb(optr,i) = Porb(iptr,i);
	    AXorb(optr,i)= AXorb(iptr,i);
	    AYorb(optr,i)= AYorb(iptr,i);
	    AYorb(optr,i)= AYorb(iptr,i);
#endif
        }
}

/*------------------------------------------------------------------------------
 * LIST_ORBIT:  simply list the orbits positions and velocities
 *------------------------------------------------------------------------------
 */
 
void list_orbit (orbitptr optr, double tstart, double tend, int n, string f)
{
    int i, kount;
    char fmt7[256];
    char fmt11[256];

    sprintf(fmt7, "%%d %s %s %s %s %s %s %s  %s %s %s %s\n",f,f,f,f,f,f,f, f,f,f,f);
    sprintf(fmt11,"%%d %s %s %s %s %s %s %s\n",             f,f,f,f,f,f,f);
        
    dprintf (0,"Total number of steps = %d\n",Nsteps(optr));
    dprintf (0,"Mass = %f \n",Masso(optr));
    dprintf (0,"Key = %d \n",Key(optr));
    dprintf (0,"Integrals of Motion (IOM) = ");
    for (i=0; i<Ndim(optr); i++)
        dprintf (0," %f ",*(IOM(optr)+i));
#ifdef PJT312
    dprintf (0," ( ",*(IOM(optr)+i));
#else
    dprintf (0," ( ");
#endif
    for (i=0; i<Ndim(optr); i++)
        dprintf (0," %f ",*(IOMERR(optr)+i));
    dprintf (0," )\n");
    dprintf (0,"Potential: Name: %s Pars: %s File: %s\n",
        PotName(optr), PotPars(optr), PotFile(optr));
    kount = 0;
    for (i=0; i<Nsteps(optr); i++) {
        if ((tstart<Torb(optr,i)) && (Torb(optr,i)<tend)) {
            if (kount++ == 0)
#ifdef ORBIT_PHI
                printf (fmt7,
			i,Torb(optr,i),Xorb(optr,i),Yorb(optr,i),Zorb(optr,i),
			Uorb(optr,i),Vorb(optr,i),Worb(optr,i));
#else	      
                printf (fmt11,
			i,Torb(optr,i),Xorb(optr,i),Yorb(optr,i),Zorb(optr,i),
			Uorb(optr,i),Vorb(optr,i),Worb(optr,i),
			Porb(optr,i),
			AXorb(optr,i),AYorb(optr,i),AZorb(optr,i));
#endif	    
            if (kount==n)
                kount=0;
        }
    }
}


#ifdef TESTBED
#include <getparam.h>

string defv[] = {
        "name=???\n    Filename",
        "mode=w\n      R or W",
        "VERSION=1.3\n 10-dec-2019 pjt",
        NULL,
};

string usage = "orbit testbed";

void ini_orbit(orbitptr);

nemo_main()
{
    orbitptr    optr=NULL;
    stream      instr,outstr;
    string      fname=getparam("name");
    
    if (strcmp(getparam("mode"),"w")==0) {              /* write test */
        printf ("WRITING test\n");
        outstr = stropen (fname,"w");
        allocate_orbit(&optr,3,4);      /* 3-D  4 steps */
        ini_orbit (optr);
        write_orbit (outstr,optr);
        strclose(outstr);
    } else {
        printf ("READING test\n");
        instr = stropen (fname,"r");
        while (read_orbit(instr,&optr)) {
                list_orbit(optr,-10.0,10.0,1,"%g");
        }
        strclose(instr);
    }
}

void ini_orbit (optr)
orbitptr optr;
{
    int i;
        
    Masso(optr)  = 1.0;
    Key(optr) = 1;
    I1(optr) = -1.0;
    I2(optr) = 0.0;
    I3(optr) = 0.0;

    PotName(optr) = "";
    PotPars(optr) = "";
    PotFile(optr) = "";

    for (i=0; i<Nsteps(optr); i++) {
        Torb(optr,i) = 0.1*i;
        Xorb(optr,i) = 0.1*i+1.0;
        Yorb(optr,i) = 0.1*i+2.0;
        Zorb(optr,i) = 0.1*i+3.0;
        Uorb(optr,i) = -Xorb(optr,i);
        Vorb(optr,i) = -Yorb(optr,i);
        Worb(optr,i) = -Zorb(optr,i);
    }
}
#endif
