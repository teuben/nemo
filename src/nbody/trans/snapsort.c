/*
 * SNAPSORT: sort particles according to a user-specified ranking.
 *
 *    22-jan-89	  V1.1   some changes	JEB
 *     9-dec-90   V1.2   helpvec	PJT
 *    16-jun-92   V1.3   usage; fixed bug with sandwiched history   PJT
 *    21-dec-92   V1.4   selectable sort routine                    PJT
 *    10-jun-95       a  declaration fix (linux) n                  pjt
 *    31-dec-02   V1.5   gcc3/SINGLEPREC                            pjt
 *     1-nov-07       a  bug when Aux is present                    pjt
 *    29-feb-08          fix a memory leak on btab                  jcl
 *    19-Jun-09          fix a bug when Aux is present              jcl
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>

#include <snapshot/body.h>
#include <snapshot/snapshot.h>
#include <snapshot/get_snap.c>
#include <snapshot/put_snap.c>
#include <bodytransc.h>

string defv[] = {	
    "in=???\n		Input file name (snapshot)",
    "out=???\n		Output file name (snapshot)",
    "rank=etot\n	Value used in ranking particles",
    "times=all\n        Range of times to process ",
    "sort=qsort\n       Sort mode {qsort;...}",
    "VERSION=1.5c\n     19-jun-09 PJT ",
    NULL,
};

string usage="sort particles according to a user-specified ranking";

string cvsid="$Id$";

/* #define FLOGGER 1       /* merge in the cute flogger test routines */

void snapsort(Body *, int , real , bool, bool, rproc_body, iproc);

nemo_main()
{
    stream instr, outstr;
    string times;
    rproc_body rank;
    iproc mysort, getsort();
    Body *btab = NULL;
    int nbody, bits;
    real tsnap;
    bool Qaux;
    

    instr = stropen(getparam("in"), "r");
    get_history(instr);
    outstr = stropen(getparam("out"), "w");
    put_history(outstr);
    rank = btrtrans(getparam("rank"));
    times = getparam("times");
    mysort = getsort(getparam("sort"));
    do {
        get_history(instr);  /*  get history that is not written */
	get_snap_by_t(instr, &btab, &nbody, &tsnap, &bits, times);
	if (bits & PhaseSpaceBit) {
	  snapsort(btab, nbody, tsnap, bits&AuxBit, bits&KeyBit, rank, mysort);
	    bits |= KeyBit;
	    put_snap(outstr, &btab, &nbody, &tsnap, &bits);
	}
	if (btab) free((Body *) btab);
	btab=NULL;	/* 'free' the snapshot */
    } while (bits != 0);
    strclose(outstr);
}


int rank_aux(Body *a, Body *b)
{
    return (Aux(a) < Aux(b) ? -1 : Aux(a) > Aux(b) ? 1 : 0);
}

void snapsort(
	 Body *btab,
	 int nbody,
	 real tsnap,
	 bool Qaux,
	 bool Qkey,
	 rproc_body rank,
	 iproc mysort)
{
    int i;
    Body *b;
    real *aux;

    if (Qaux) {  /* make backup copy of Aux */
      aux = (real *) allocate(nbody*sizeof(real));
      for (i = 0, b = btab; i < nbody; i++, b++)
	aux[i] = Aux(b);
    }
    if (!Qkey) { /* initialize Key's if they didn't exist */
      for (i = 0, b = btab; i < nbody; i++, b++)
	Key(b) = i;
    }
    

    for (i = 0, b = btab; i < nbody; i++, b++) {
	Aux(b) = (rank)(b, tsnap, i);
    }
    (mysort)(btab, nbody, sizeof(Body), rank_aux);

    if (Qaux) {   /* stuff it back */
      for (i = 0, b = btab; i < nbody; i++, b++)
	Aux(b) = aux[Key(b)];

      free(aux);
    }
}



/*
 *  Tabulate the valid sort names, plus their associated external
 *  routines. The accompanying getsort() routine returns the 
 *  appropriate sort routine
 */

typedef struct sortmode {
    string name;
    iproc   fie;
} sortmode;

#define SortName(x)  x->name
#define SortProc(x)  x->fie


/* List of externally available sort routines */

/* extern int qsort();         /* Standard Unix : stdlib.h */
/* or:  void qsort(void *base, size_t nmemb, size_t size,
 *                 int(*compar)(const void *, const void *));
 */
                  
extern int bubble_sort();          /* Flogger library routines */
extern int heap_sort();
extern int insertion_sort();
extern int merge_sort();
extern int quick_sort();
extern int shell_sort();


local sortmode smode[] = {
#ifdef FLOGGER
    "bubble",   bubble_sort,        /* cute flogger routines */
    "heap",     heap_sort,
    "insert",   insertion_sort,
    "merge",    merge_sort,
    "quick",    quick_sort,
    "shell",    shell_sort,
#endif
    "qsort",    (iproc) qsort,      /* standard Unix qsort() */
    NULL, NULL,
};

iproc getsort(string name)
{
    sortmode *s;

    for (s=smode; SortName(s); s++)  {
        dprintf(1,"GETSORT: Trying %s\n",SortName(s));
        if (streq(SortName(s),name))
            return SortProc(s);
    }
    error("No valid sortname");
    return NULL;     /* better not get here ... */
}
