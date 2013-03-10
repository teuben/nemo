/*
 * RAN3: loosely based on ran3() from NumRec
 *
 *	17-sep-95: addition safeguard against out of bounds
 *		   when initialized with large idum's (idum > MBIG)
 *       9-may-97: option to view long as int for 64bit machines to keep
 *                 it portable
 *      13-may-05: made int the default, which works fine as long
 *                 as sizeof(int) = 4; 
 *                 Also fixed the prototype for ran3()
 *
 *  Although based on his original version, Knuth has since then published
 *  a much improved version of ran3() in
 *  "The Stanford GraphBase," Addison-Wesley, 1993 (pp 216-221).
 *
 */


#include <stdinc.h>

#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

#define REDIAL      /* Set this if you want to check for inbounds & redial */


#if 1
typedef int mylong;
#else
typedef long mylong;  /* this is the original definition, but may fail on IA64 */
#endif

real ran3(int *idum)
{
    /* -- save inbetween calls */
    local int inext,inextp;
    local mylong ma[56];		/* never modify this !! (see Knuth) */
    local int iff=0;
    /* -- scratch variables */
    mylong mj,mk;
    int i,ii,k;

    if (*idum < 0 || iff == 0) {		/* initialize */
        if (sizeof(mylong) != 4) 
	    warning("ran3: may not be portable: sizeof(mylong)=%d",sizeof(mylong));
        iff=1;
	mj=MSEED-(*idum < 0 ? -*idum : *idum);
#ifdef REDIAL
        while(mj < MZ)mj += MBIG;    
#endif		
	mj %= MBIG;
	ma[55]=mj;
	mk=1;
	for (i=1;i<=54;i++) {			/* initialize rest of table */
	    ii=(21*i) % 55;
	    ma[ii]=mk;
	    mk=mj-mk;
#ifdef REDIAL
	    while (mk < MZ) mk += MBIG;
#else			
	    if (mk < MZ) mk += MBIG;
#endif			
	    mj=ma[ii];
	}
	for (k=1;k<=4;k++)
	    for (i=1;i<=55;i++) {
		ma[i] -= ma[1+(i+30) % 55];
		if (ma[i] < MZ) ma[i] += MBIG;
	    }
	inext=0;
	inextp=31;				/* special 'Knuth' number */
	*idum=1;
    }
    if (++inext == 56) inext=1;
    if (++inextp == 56) inextp=1;
    mj=ma[inext]-ma[inextp];
#ifdef REDIAL    
    while (mj < MZ) mj += MBIG;
#else    
    if (mj < MZ) mj += MBIG;
#endif    
    ma[inext]=mj;
    return mj*FAC;
}

/* (C) Copr. 1986-92 Numerical Recipes Software ?421.1-9. */
