/*
 * TABPAIRS:  star encounter histories; see also man page
 *           
 *      This program reads a grep/awk trimmed table from a starlab run,
 *      and performs various statistics on this table.
 *
 *      
 *
 *	21-jun-00	PJT @ AMNH - figuring out pieces of star encounter histo
 *       1-mar-03       fixed file_lines
 *       1-jan-04       get_line changed interface
 *
 */

#include <stdinc.h>	
#include <getparam.h>
#include <moment.h>
#include <extstring.h>

#define MLINELEN    512
#define MAXLIST       8
#define MAXPAIRS 150000
#define MAXSTARS  10000

string defv[] = {                /* DEFAULT INPUT PARAMETERS */
    "in=???\n            Input file name (table)",
    "star=\n             Star(s) to follow",
    "count=f\n           Count occurance of all stars?",
    "delete=0\n		 Iteration count to delete obvious combo/split pairs?",
    "nmax=10000\n        maximum number of data to be read if pipe",
    "VERSION=1.0b\n	 1-jan-04 PJT",
    NULL
};

string usage = "table histories";


local string input;				/* filename */
local stream instr;				/* input file */

local int star;
local int nmax;
local bool Qcount;
local int iter_del;

void setparams(void), isort (int *x , int n), read_data(void),
     stat_data(void), star_data(int), count_data(void), delete_data(void);


extern string *burststring(string,string);


typedef struct {   /* some convenient structure to hold all pair (encounters) in memory */
  int k;                  /* k=1 combo,  k=-1 split */
  int n;                  /* length of nlist */
  int nlist[MAXLIST];     /* star numbers (0 based) */
  real tsnap;             /* time */
  real e;                 /* energy of "binary" */
  real p;                 /* period, if e < 0 */
} Pair;


Pair pair[MAXPAIRS];
int npairs;

int combo_cnt[MAXSTARS], split_cnt[MAXSTARS];


void nemo_main()
{
    setparams();

    read_data();
    delete_data();
    if (star > 0) {
        star_data(star);
    } else if (Qcount) {
        count_data();
    } else
        stat_data();
}

void setparams(void)
{
    input = getparam("in");             /* input table file */
    instr = stropen (input,"r");

    if (hasvalue("star"))
        star = getiparam("star");
    else
        star = -1;

    Qcount = getbparam("count");
    iter_del = getiparam("delete");
    
    nmax = nemo_file_lines(input,getiparam("nmax"));
    dprintf(0,"Allocated %d lines for table\n",nmax);

}

/*
 * Need to read lines like this:
 *
 *  1 7.83933 -0.00583491 0.996812 5574 5995
 * -1 7.8396 -0.00583602 0.996528 (5995,5574)
 *
 */
 
void read_data(void)
{
    int i, k, nwords, n, n1, n2, count=0;
    char line[MLINELEN];
    real tsnap, e, p;
    string *sp, *nlist1, *nlist2;
    int nlist[MAXLIST];
		
    for(;;) {
        if (get_line(instr,line) < 0) {
            dprintf(1,"## Read %d lines\n",count);
            npairs = count;
            break;
        }
        count++;

        sp = burststring(line," ");
        nwords = xstrlen(sp,sizeof(string))-1;
        k = atoi(sp[0]);
        tsnap = atof(sp[1]);
        e = atof(sp[2]);
        p = atof(sp[3]);
        nlist1 = burststring(sp[4],",()");
        n1 = xstrlen(nlist1,sizeof(string))-1;
        for (i=0, n=0; i<n1; i++)
            nlist[n++] = atoi(nlist1[i]);
        if (k > 0) {        /* creating a binary */
            nlist2 = burststring(sp[5],",()");
            n2 = xstrlen(nlist2,sizeof(string))-1;
        } else
            n2 = 0;
        for (i=0; i<n2; i++)
            nlist[n++] = atoi(nlist2[i]);

        dprintf(1,"%d : %d %d :",count,n1,n2);
        for (i=0; i<n; i++)
            dprintf(1," %d",nlist[i]);
        dprintf(1,"\n");

        isort(nlist,n);
#if 0        
        printf("%s %s %s %s",sp[0],sp[1],sp[2],sp[3]);
        for (i=0; i<n; i++)
            printf(" %d",nlist[i]);
        printf("\n");
#endif

        pair[count-1].k = k;
        pair[count-1].n = n;
        pair[count-1].tsnap = tsnap;
        pair[count-1].e = e;
        pair[count-1].p = p;
        for (i=0; i<n; i++)
            pair[count-1].nlist[i] = nlist[i];
	for (i=n; i<MAXLIST; i++)
            pair[count-1].nlist[i] = -1;
    }
}


void old_stat_data(void)
{
    int i,j, n, notsame;
    
    
    for (i=1; i<npairs; i++) {
        if (pair[i].n == pair[i-1].n) {
            n = pair[i].n;
            notsame = 0;
            for (j=0; j<n; j++)
                if (pair[i].nlist[j] != pair[i-1].nlist[j]) notsame++;
            if (notsame == 0) {     /* ok, they are same ; don't print */
                i++;
            }
        } else
            notsame = 1;

        if (notsame) {
            printf("%d %g %g %g", i+1, pair[i].tsnap, pair[i].e,pair[i].p);
            for (j=0; j<pair[i].n; j++)
                printf(" %d",pair[i].nlist[j]);
            printf("\n");
	}            

    }
}

void stat_data(void)
{
    int i,j;
    
    
    for (i=0; i<npairs; i++) {
        if (pair[i].n > 0) {
            printf("%d %g %g %g", i+1, pair[i].tsnap, pair[i].e,pair[i].p);
            for (j=0; j<pair[i].n; j++)
                printf(" %d",pair[i].nlist[j]);
            printf("\n");
	}            
    }
}

void star_data(int star)
{
    int n, i, j, gotit;
    
    for (i=0; i<npairs; i++) {
        n = pair[i].n;
        if (n==0) continue;
        gotit = 0;
        for (j=0; j<n; j++) {
            if (pair[i].nlist[j] == star) {
                gotit=1;
                break;
            }
        }
        if (gotit) {
            printf("%d %d %g %g %g",i+1,
			pair[i].k, pair[i].tsnap, pair[i].e,pair[i].p);
            for (j=0; j<pair[i].n; j++)
                printf(" %d",pair[i].nlist[j]);
            printf("\n");
	}            
    }

}


/*
 *  This should remove (set nlist length (pair.n) to zero) all entries
 *  which are obvious pairs of combo/split with the same interaction list
 */
void delete_data(void)
{
    int i,j, n, notsame, ndel=0;
    int sep=1;
    bool Qdel = (iter_del > 0);

    for (sep=1; sep <= iter_del; sep++) {

        for (i=sep; i<npairs; i++) {
            n = pair[i].n;
            if (n>0  && n == pair[i-sep].n) {
                notsame = 0;
                for (j=0; j<n; j++) 
                    if (pair[i].nlist[j] != pair[i-sep].nlist[j]) notsame++;
                if (notsame == 0) {     /* ok, they are same ; delete */
                    if (Qdel) {
                        pair[i].n = pair[i-sep].n = 0;
                        ndel++;
                    }
                }
            } 
    	}
    	if(Qdel) dprintf(0,"Deleted %d obvious combo/split pairs\n",ndel);
    }
}

void count_data(void)
{
    int i,j,k,n, mstar=0;

    dprintf(1,"Counting.....\n");

    for (i=0; i<MAXSTARS; i++)
        combo_cnt[i] = split_cnt[i] = 0;
    
    for (i=0; i<npairs; i++) {
        n = pair[i].n;
        for (j=0; j<n; j++) {
            k = pair[i].nlist[j];
            if (k>MAXSTARS) error("star number %d too high",k);
            mstar = MAX(mstar,k);
            if (pair[i].k > 0)
                combo_cnt[k]++;
            else
                split_cnt[k]++;
        }
    }
    dprintf(1,"Max star # found: %d\n",mstar);

    for (i=0; i<=mstar; i++) {
        if (combo_cnt[i]>0 || split_cnt[i]>0)
            printf("%d %d %d\n",i,combo_cnt[i],split_cnt[i]);
    }


}




void isort (int *x , int n)
{
    int    gap, i, j, temp;

    if (n==2) {
    	if (x[0] > x[1]) {
            temp = x[0];
            x[0] = x[1];
            x[1] = temp;
    	}
    	return;
    }

    for (gap=n/2; gap>0; gap /= 2)
        for (i=gap; i<n; i++)
            for (j=i-gap; j>=0; j -= gap) {
                if (x[j] <= x[j+gap])
                    break;
                temp = x[j];
                x[j] = x[j+gap];
                x[j+gap] = temp;

            }
}

