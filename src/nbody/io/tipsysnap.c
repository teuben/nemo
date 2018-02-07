/* 
 * TIPSYSNAP: convert ascii tipsy to snapshot(5NEMO)
 *
 * $Header$
 * $Log$
 * Revision 1.4  2003/10/25 15:47:42  pteuben
 * moved declration
 *
 * Revision 1.3  2001/08/28 21:37:06  pteuben
 * fix reading binary data, got rid of hardcoded TIPSY_NEEDPAD
 * and warn when header size is not 28.
 *
 * Revision 1.2  2001/04/02 05:18:51  teuben
 * more fixing of the linux install
 *
 * Revision 1.1.1.1  2000/08/19 03:56:44  teuben
 * import NEMO V3
 *
 * Revision 1.1  93/09/14  09:51:22  trq
 * Initial revision
 * 
 * Revision 2.2  91/06/02  12:21:32  trq
 * Added star stuff from Neal
 * 
 * Revision 2.1  1991/04/18  18:11:35  trq
 * Added checks for malloc(0).
 *
 *  2-sep-94	For NEMO conversion to snapshot
 * 18-aug-00    also allow binary files..., with padding and optional swap
 *  1-apr-01    compiler warnings - pjt
 * 26-aug-01    fix problems with mass99 data that have nsph=ndark=0
 */
 
#include <stdinc.h>
#include <getparam.h>
#include <math.h>
#include <stdlib.h>
#include <vectmath.h>		/* otherwise NDIM undefined */
#include <filestruct.h>
#include <history.h>

#include <snapshot/snapshot.h>	
#include <snapshot/body.h>
#include <snapshot/put_snap.c>
#include "tipsydefs.h"              /* header is in here */

string defv[] = {
    "in=???\n                   Input (tipsy) ascii file",
    "out=???\n                  Output snapshot file",
    "options=gas,dark,star\n    Output which particles?",
    "mode=ascii\n		Input mode (ascii, binary)",
    "swap=f\n                   Swap bytes?",
    "offset=0\n                 Offset data from header?",
    "VERSION=2.1b\n             4-feb-2018 pjt",
    NULL,
};

string usage="convert tipsy ascii/binary file to snapshot";

void nemo_main()
{
    int ndim ;
    int nbodies ;
    int ngas ;
    int ndark ;
    int nstar ;
    int last_nstar ;
    int count ;
    float dummy ;
    float *tform  = NULL;
    float *tf ;
    float last_tform ;
    struct  gas_particle *gp, *lastgp;
    struct dark_particle *dp, *lastdp ;
    struct star_particle *sp, *lastsp, *last_old_sp ;

    /* NEMO */
    bool Qgas, Qdark, Qstar, Qascii, Qbinary, Qswap;
    string options, mode;
    int i, j, n, nbody, bits, offset;
    real tsnap;
    Body *btab = NULL, *bp;
    stream instr, outstr;
    
    instr = stropen(getparam("in"),"r");
    outstr = stropen(getparam("out"),"w");
    put_history(outstr);
    
    options = getparam("options");
    Qgas  = scanopt(options,"gas");
    Qdark = scanopt(options,"dark");
    Qstar = scanopt(options,"star");

    mode = getparam("mode");
    Qascii   = scanopt(mode,"ascii");
    Qbinary  = scanopt(mode,"binary");

    Qswap = getbparam("swap");
    
    offset = getiparam("offset");

    last_nstar = 0 ;
    last_tform = -1.e30 ;
    if (Qascii) for(;;) {
	count=fscanf(instr, "%d%*[, \t\n]%d%*[, \t\n]%d",
		    &header.nbodies, &header.nsph, &header.nstar) ;
	if ( (count == EOF) || (count==0) )  break ;

	fscanf(instr,"%d",&header.ndim) ;
	fscanf(instr,"%lf",&header.time) ;
	ndim=header.ndim;
	nbodies=header.nbodies;
	ngas=header.nsph;
	nstar = header.nstar ;
	ndark = header.ndark = nbodies - nstar - ngas ;

        nbody = 0;
        if (Qgas) nbody += ngas;
        if (Qdark)nbody += ndark;
        if (Qstar)nbody += nstar;
        tsnap = header.time;
	dprintf(0,"read time %lf input (%d,%d,%d) output %d\n",
	    header.time, ngas, ndark, nstar, nbody) ;
        if (nbody==0) error("Nothing to output");

	if(gas_particles != NULL) free(gas_particles);
	if(ngas != 0)
	  gas_particles = (struct gas_particle *) 
	            allocate(ngas*sizeof(*gas_particles));
	else
	  gas_particles = NULL;
	  
	if(dark_particles != NULL) free(dark_particles);
	if(ndark != 0)
	  dark_particles = (struct dark_particle *) 
	            allocate(ndark*sizeof(*dark_particles));
	else
	  dark_particles = NULL;
	  
	if(star_particles != NULL) free(star_particles);
	if(nstar != 0) {
	    star_particles = (struct star_particle *)
	        allocate(nstar*sizeof(*star_particles));
	    if(tform == NULL)
	      tform = (float *)allocate(nstar*sizeof(*tform));
	    else if(nstar > last_nstar)
	      tform = (float *)reallocate(tform,nstar*sizeof(*tform));
	} else
	  star_particles = NULL;
	  
        if (btab != NULL) free(btab);
        btab = (Body *) allocate(nbody*sizeof(Body));
	 
	lastgp = gas_particles + ngas ;
	lastdp = dark_particles + ndark ;
	lastsp = star_particles + nstar ;
	last_old_sp = star_particles + last_nstar ;

	for(gp=gas_particles; gp < lastgp ; gp++)
	    fscanf(instr,"%f%*[, \t\n]",&gp->mass);
	for(dp=dark_particles; dp < lastdp;  dp++)
	    fscanf(instr,"%f%*[, \t\n]",&dp->mass);
	for(sp=star_particles; sp < lastsp; sp++)
	    fscanf(instr,"%f%*[, \t\n]",&sp->mass);
	for(gp=gas_particles; gp < lastgp ; gp++)
	    fscanf(instr,"%f%*[, \t\n]",&gp->pos[0]);
	for(dp=dark_particles; dp < lastdp ; dp++) 
	    fscanf(instr,"%f%*[, \t\n]",&dp->pos[0]);
	for(sp=star_particles; sp < lastsp ; sp++) 
	    fscanf(instr,"%f%*[, \t\n]",&sp->pos[0]);
	for(gp=gas_particles; gp < lastgp ; gp++) 
	    fscanf(instr,"%f%*[, \t\n]",&gp->pos[1]);
	for(dp=dark_particles; dp < lastdp ; dp++) 
	    fscanf(instr,"%f%*[, \t\n]",&dp->pos[1]);
	for(sp=star_particles; sp < lastsp ; sp++) 
	    fscanf(instr,"%f%*[, \t\n]",&sp->pos[1]);
	if (ndim == 3){
	    for(gp=gas_particles; gp < lastgp ; gp++)
		fscanf(instr,"%f%*[, \t\n]",&gp->pos[2]);
	    for(dp=dark_particles; dp < lastdp ; dp++) 
		fscanf(instr,"%f%*[, \t\n]",&dp->pos[2]);
	    for(sp=star_particles; sp < lastsp ; sp++)
		fscanf(instr,"%f%*[, \t\n]",&sp->pos[2]);
	}
	for(gp=gas_particles; gp < lastgp ; gp++)
	    fscanf(instr,"%f%*[, \t\n]",&gp->vel[0]);
	for(dp=dark_particles; dp < lastdp ; dp++) 
		fscanf(instr,"%f%*[, \t\n]",&dp->vel[0]);
	for(sp=star_particles; sp < lastsp ; sp++) 
	    fscanf(instr,"%f%*[, \t\n]",&sp->vel[0]);
	for(gp=gas_particles; gp < lastgp ; gp++)
	    fscanf(instr,"%f%*[, \t\n]",&gp->vel[1]);
	for(dp=dark_particles; dp < lastdp ; dp++) 
	    fscanf(instr,"%f%*[, \t\n]",&dp->vel[1]);
	for(sp=star_particles; sp < lastsp ; sp++) 
	    fscanf(instr,"%f%*[, \t\n]",&sp->vel[1]);
	if (ndim == 3){
	    for(gp=gas_particles; gp < lastgp ; gp++) 
		fscanf(instr,"%f%*[, \t\n]",&gp->vel[2]);
	    for(dp=dark_particles; dp < lastdp ; dp++) {
		count = fscanf(instr,"%f%*[, \t\n]",&dp->vel[2]);
		if ( (count == EOF) || (count==0) )
		    error("short read on vz");
	    }
	    for(sp=star_particles; sp < lastsp ; sp++)
		fscanf(instr,"%f%*[, \t\n]",&sp->vel[2]);
	}
	for(dp=dark_particles; dp < lastdp ; dp++) {
	    count = fscanf(instr,"%f%*[, \t\n]",&dp->eps);
	    if ( (count == EOF) || (count==0) ) 
		error("short read on epsilons");
	}
	for(sp=star_particles; sp < lastsp ; sp++)
	    fscanf(instr,"%f%*[, \t\n]",&sp->eps);
	for(gp=gas_particles; gp < lastgp ; gp++)
	    fscanf(instr,"%f%*[, \t\n]",&gp->rho);
	for(gp=gas_particles; gp < lastgp ; gp++)
	    fscanf(instr,"%f%*[, \t\n]",&gp->temp);
	for(gp=gas_particles; gp < lastgp ; gp++)
	    fscanf(instr,"%f%*[, \t\n]",&gp->hsmooth);
	for(gp=gas_particles; gp < lastgp ; gp++)
	    fscanf(instr,"%f%*[, \t\n]",&gp->metals);
	for(sp=star_particles; sp < lastsp ; sp++) 
	    fscanf(instr,"%f%*[, \t\n]",&sp->metals);
	for(sp=star_particles, tf = tform ; sp < last_old_sp ; sp++, tf++) {
	    fscanf(instr,"%f%*[, \t\n]",&dummy);
	    sp->tform = *tf ;
	}
	for(sp=star_particles + last_nstar, tf = tform + last_nstar ;
		sp < lastsp ; sp++,tf++) {
	    fscanf(instr,"%f%*[, \t\n]",tf);
	    if(*tf < last_tform){
	        warning("error in star formation time: %g < %g", *tf,last_tform); 
		return ;
	    }
	    sp->tform = *tf ;
	    last_tform = *tf ;
	}
	last_nstar = nstar ;
	for(gp=gas_particles; gp < lastgp ; gp++){
	    count = fscanf(instr,"%f%*[, \t\n]",&gp->phi);
	    if ( (count == EOF) || (count==0) )
	      return ;
	}
	for(dp=dark_particles; dp < lastdp;  dp++){
	    count = fscanf(instr,"%f%*[, \t\n]",&dp->phi);
	    if ( (count == EOF) || (count==0) )
		error("short read on potentials");
	}
	for(sp=star_particles; sp < lastsp; sp++){
	    count = fscanf(instr,"%f%*[, \t\n]",&sp->phi);
	    if ( (count == EOF) || (count==0) )
	      return ;
	}
	bp=btab;
	if (Qgas)	    
	    for (i=0; i<ngas; i++) {
	        Mass(bp) = gas_particles[i].mass;
	        for (j=0; j<ndim; j++) {
	            Pos(bp)[j] = gas_particles[i].pos[j];
	            Vel(bp)[j] = gas_particles[i].vel[j];
	        }
	        bp++;
	    }
        if (Qdark)
	    for (i=0; i<ndark; i++) {
	        Mass(bp) = dark_particles[i].mass;
	        for (j=0; j<ndim; j++) {
	            Pos(bp)[j] = dark_particles[i].pos[j];
	            Vel(bp)[j] = dark_particles[i].vel[j];
	        }
	        bp++;
	    }
        if (Qstar)
	    for (i=0; i<nstar; i++) {
	        Mass(bp) = star_particles[i].mass;
	        for (j=0; j<ndim; j++) {
	            Pos(bp)[j] = star_particles[i].pos[j];
	            Vel(bp)[j] = star_particles[i].vel[j];
	        }
	        bp++;
	    }
    
        bits = (TimeBit | MassBit | PhaseSpaceBit);	    
        put_snap(outstr, &btab, &nbody, &tsnap, &bits);
    }


    if (Qbinary) for(;;) {
        n = fread((char *)&header,sizeof(header),1,instr) ;
        if (n <= 0) break;
	if (sizeof(header) != 28)
	  warning("TIPSY header = %d bytes, not 28. You may need offset=",sizeof(header));
#ifdef TIPSY_NEEDPAD
	dprintf(1,"TIPSY_NEEDPAD was defined; header=%d bytes\n",sizeof(header));
#else
	dprintf(1,"TIPSY_NEEDPAD was not defined; header=%d bytes\n",sizeof(header));
#endif	
        if (Qswap) {
            bswap((void *)&header.time,    sizeof(double), 1);
            bswap((void *)&header.nbodies, sizeof(int),    1);
            bswap((void *)&header.ndim,    sizeof(int),    1);
            bswap((void *)&header.nsph,    sizeof(int),    1);
            bswap((void *)&header.ndark,   sizeof(int),    1);
            bswap((void *)&header.nstar,   sizeof(int),    1);
        }

	ndim=header.ndim;

	nbodies=header.nbodies;
	ngas=header.nsph;
	nstar = header.nstar ;
	ndark = header.ndark = nbodies - nstar - ngas ;

        switch(header.ndim) {
            case 3:
                    dprintf(0,"ndim=%d\n",header.ndim);
                    break;
            default:
                    error("strange header.ndim=%d, need swap?",header.ndim);
                    break;
        }

        dprintf(0,"time=%g N=%d (%d,%d,%d) for gas,dark,star\n",
                header.time, header.nbodies,
                header.nsph, header.ndark, header.nstar);

	if (offset) {
	  n = fseek(instr, offset, SEEK_CUR);
	  warning("Stupid TIPSY, offsetting data by %d bytes....: %d",n);
	}

        if (header.nsph > 0) {
            gp = allocate(header.nsph*sizeof(*gp));
            n = fread((char *)gp, sizeof(*gp), header.nsph, instr);
            if (n <= 0) error("Problem reading gas data");            
	    dprintf(1,"Found %d gas particles\n",n);

	    if (Qgas) {
	      nbody = ngas;
	      if (btab != NULL) free(btab);
	      btab = (Body *) allocate(nbody*sizeof(Body));
	      bp=btab;
	      for (i=0; i<nbody; i++) {
	        Mass(bp) = gp[i].mass;
	        for (j=0; j<ndim; j++) {
		  Pos(bp)[j] = gp[i].pos[j];
		  Vel(bp)[j] = gp[i].vel[j];
	        }
	        bp++;
	      }
	      bits = (TimeBit | MassBit | PhaseSpaceBit);	    
	      put_snap(outstr, &btab, &nbody, &tsnap, &bits);

	    }
        }

        if (header.ndark > 0) {
            dp = allocate(header.ndark*sizeof(*dp));
            n = fread((char *)dp, sizeof(*dp), header.ndark, instr);
            if (n <= 0) error("Problem reading dark data");            
	    dprintf(1,"Found %d dark matter particles\n",n);

	    if (Qdark) {
	      nbody = ndark;
	      if (btab != NULL) free(btab);
	      btab = (Body *) allocate(nbody*sizeof(Body));
	      bp=btab;
	      for (i=0; i<nbody; i++) {
	        Mass(bp) = dp[i].mass;
	        for (j=0; j<ndim; j++) {
		  Pos(bp)[j] = dp[i].pos[j];
		  Vel(bp)[j] = dp[i].vel[j];
	        }
	        bp++;
	      }
	      bits = (TimeBit | MassBit | PhaseSpaceBit);	    
	      put_snap(outstr, &btab, &nbody, &tsnap, &bits);
	    }
        }

        if (header.nstar > 0) {
            sp = allocate(header.nstar*sizeof(*sp));
            n = fread((char *)sp, sizeof(*sp), header.nstar, instr);
            if (n <= 0) error("Problem reading star data");            
	    dprintf(1,"Found %d star particles\n",n);

	    if (Qstar) {
	      nbody = nstar;
	      if (btab != NULL) free(btab);
	      btab = (Body *) allocate(nbody*sizeof(Body));
	      bp=btab;
	      for (i=0; i<nbody; i++) {
	        Mass(bp) = sp[i].mass;
	        for (j=0; j<ndim; j++) {
		  Pos(bp)[j] = sp[i].pos[j];
		  Vel(bp)[j] = sp[i].vel[j];
	        }
	        bp++;
	      }
	      bits = (TimeBit | MassBit | PhaseSpaceBit);	    
	      put_snap(outstr, &btab, &nbody, &tsnap, &bits);

	    }
        }
    }

    if(tform) free(tform);
    strclose(instr);
    strclose(outstr);
}
