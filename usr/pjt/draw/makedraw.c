/*
  MakeDraw:   make a draw, starting from a ranked list of teams
 
  The list of teams should be sorted, such that #1 is at the top,
  ranking is done by a parameters that controls from which location
  onwards players are considered all equal.
  The normal way is to assign a 1st and 2nd seed, the #3/4 seeds are
  considered equal, and 5-8 can be considered shadow seeds. Players
  9-16 can again be clumped together, but considered equal.
  
 
 *
 *	0.1	12-oct-1994	initial writeup for DC Open 94 - peter teuben
 *	0.2	 2-oct-1999	cleanup for DC Open 99
 *      0.3     26-oct-1999     also allow finishing random draw (full=t)
 *      0.4     31-oct-1999     more work on the 9...16 places
 *      0.5      2-nov-1999     allow randomize on some levels
 *	0.6      4-nov-1999     optional add match IDs to matches per round
 *      0.7     25-oct-2001     matchid -1 won't write that column
 *      0.8      7-mar-2002     IBF/USAB rules
 *      0.9      5-jun-2004     minor improvements at US Nationals 2004 (Shreveport)
 *
 * Todo: 
        implement shadow and normal seeds properly, use an array of nrandom= ??
	(wo)men doubles need alphabetized order of players

   nrandom:  should this be done in blocksteps of powers of 2 ?

   Buga
   June 2004:   there are some bugs in this draw. For a draw of
                size N, the N/2+1..N slots are ok, but via
                recursion earlier draws will be messed up
                Due to the way seeding goes, e.g. 5/8 seeds
                are randomized, it will not make a difference
                in the draw.
		The problem is not in the setup_byes() routine,
		that one appears to be correct.
 */


#include <nemo.h>
#include <extstring.h>
#include <ctype.h>
#include "draw.h"

string defv[] = {
    "in=\n          Input teams for this draw",
    "out=-\n        Output draw, can be used in e.g. drawplot",
    "p=4\n          level depth ; 2**n is number of players ",
    "seed=0\n       Random seed",
    "full=f\n       Finish a full random draw (testing)",
    "auto=f\n       Automatically advance first round byes",
    "nrandom=0\n    Set slot beyond which to randomize [1..2**p]",
    "matchid=\n     Add match id ; one offset per round",
    "maxplayer=0\n  Override maxplayer for a long list (for testing)",
    "VERSION=0.9\n  7-jun-04 PJT",
    NULL,
};

string usage="Make a drop down draw from a (ranked) list of players";

extern string *burststring(string, string);
extern double  xrandom(double, double);

int     guess_p(string);
int     parse_team(string, int, int, Team *, int);
//int     rank_cmp(const Team *, Team *);
//int     slot_cmp(const Team *, const Team *);
int     rank_cmp(const void *, const void *);
int     slot_cmp(const void *, const void *);
void    make_draw(int, Team *, int, int, int);
void    make_draw_old(int, Team *, int, int, int);
void    flip(int,int *);
void    shuffle(int,int *);
void    setup_byes(int);

static int **byes;       // will be created by setup_byes()

void nemo_main()
{
    int n, i, i0, idx, p = getiparam("p");
    int nid, mid, midx, matchid[256];
    char smatchid[128];
    Team *team, *w, zero;
    stream outstr;
    int nteam, slot, maxplayer=getiparam("maxplayer");
    int seed, nrandom;
    bool Qfull = getbparam("full");
    bool Qauto = getbparam("auto");

    set_xrandom(getiparam("seed"));

    if (p==0) {    /* try and figure it out from the input list */
      p = guess_p(getparam("in"));
      if (p==0) error("Cannot guess input size, need to specify p=");
      warning("p=%d guessed %s",p,getparam("in"));
    }

    for (n=1, i=0; i<p; i++)       /* count max size of the draw */
        n *= 2;
    setup_byes(p);

    nrandom = getiparam("nrandom");
    if (hasvalue("matchid")) {
        nid = nemoinpi(getparam("matchid"),matchid,10);
        if (nid < 0) error("Matchid parsing");
    } else
        nid = 0;
    
    team = (Team *) allocate(2*n*sizeof(Team));     /* allocate twice for the whole draw !! */
    nteam = parse_team(getparam("in"),p,n, team, maxplayer);   
    dprintf(0,"Found %d players, using %d slots in the first round\n",nteam,n);
    if (nteam > n)
      error("Too many players, try setting larger p=%d",p);

    make_draw(nteam, team, p, n, nrandom);

    qsort(team,nteam,sizeof(Team),slot_cmp);
    for (i=0; i<nteam; i++) {
        dprintf(1,"Draw: %d/%d %d %d %d %s\n",i+1,n,
                team[i].id, team[i].rank, team[i].slot, team[i].name);
    }

    outstr = stropen(getparam("out"),"w");
    fprintf(outstr,"## p=%d\n",p);
    fprintf(outstr,"# First round of %d at 0 with i0=0\n",n);
    for (i=0, slot=1; i<nteam; i++, slot++) {
        while (slot < team[i].slot) {
            fprintf(outstr,"# %d/%d BYE\n",slot ,n);
            slot++;
        }
        fprintf(outstr,"%d/%d %s\n",team[i].slot ,n, team[i].name);
    }

    /* make a dummy 'zero' entry in the draw - it represents a BYE */
    zero.name[0] = 0;

    /* fill out the rest of the blank draw */
    /* can optionally add some match ID's here */

    idx = 0;                    /* index into the team's array */
    i0 = 0;                     /* offset counter for empty first round spots */
    midx = 0;                   /* index into matchid */
    smatchid[0] = 0;
    while (n > 1) {                     /* loop while still to play */
        if (nid) mid = matchid[midx++];
        if (n==2)
            fprintf(outstr,"# and the winner \n");
        else
            fprintf(outstr,"# round of %d at %d with i0=%d\n",n/2,idx+n,i0);
        for (i=0; i<n; i+=2) {                    // loop in pairs

	  if (idx == 0 && nteam<n) {              /* first round, and unfilled */

                dprintf(1,"%d: looking at %d=?%d %d=?%d\n",
                        i,
                        team[i0].slot, i+1,
                        team[i0+1].slot, i+2);
                        

                if (team[i0].slot != i+1) {                    /* upper missing */
                    dprintf(1,"Missing upper slot %d = %d\n",i+1,team[i0].slot);
                    if (team[i0].slot != i+2)  {
                        warning("two BYE matches in the same slot");
                        w = &zero;
                    } else {
                        w = &team[i0];
                        i0++;
                    }
                } else if( team[i0+1].slot != i+2)  {          /* lower missing */
                    dprintf(1,"Missing lower slot %d = %d\n",i+2,team[i0+1].slot);
                    w = &team[i0];
                    i0++;
                } else {                                        /* both present */
                    if (strlen(team[i0].name) == 0)
                        w = &team[i0+1];
                    else if (strlen(team[i0+1].name) == 0)
                        w = &team[i0];
                    else 
                        if (xrandom(-1.0,1.0) > 0)
                            w = &team[i0];
                        else
                            w = &team[i0+1];
                    i0 += 2;
                }
	     } else {  // full draw OR not the first round
                /* decide battle between idx+{i & i+1} */
	       if (Qauto) {
		 w = &zero;
	       } else {
		 if (strlen(team[idx+i].name) == 0)
                    w = &team[idx+i+1];
		 else if (strlen(team[idx+i+1].name) == 0)
                    w = &team[idx+i];
		 else 
                    if (xrandom(-1.0,1.0) > 0)
		      w = &team[idx+i];
                    else
		      w = &team[idx+i+1];
	       }
            }
            if (nid) {
	      if (mid > 0) 
		sprintf(smatchid,"[##%d] ",mid++);
	      else
		sprintf(smatchid," ");
	    }

            strcpy(team[idx+n+i/2].name, w->name);
            team[idx+n+i/2].slot = idx + n + (i + 2)/2;
            if (Qfull) {
                if(strlen(team[idx+n+i/2].name)==0)
                    fprintf(outstr,"# %d/%d BYE\n",(i+2)/2,n/2);
                else
                    fprintf(outstr,"%d/%d %s\n",(i+2)/2,n/2,team[idx+n+i/2].name);
            } else if (Qauto) {
                if(strlen(team[idx+n+i/2].name)==0)
		    fprintf(outstr,"%d/%d %s\n",(i+2)/2,n/2,smatchid);
		else
                    fprintf(outstr,"%d/%d %s\n",(i+2)/2,n/2,team[idx+n+i/2].name);
	    } else
                fprintf(outstr,"%d/%d %s\n",(i+2)/2,n/2,smatchid);
		
        }
        idx += n;
        n /= 2;
    }
    strclose(outstr);

}


#define MAXLINELEN 256

int guess_p(string fname)
{
  stream instr;
  char line[MAXLINELEN], *cp;
  int len;
  int  np, n=0, p=0;

  instr = stropen(fname,"r");  
  while (fgets(line,MAXLINELEN,instr) != NULL) {
    len = strlen(line);
    cp = &line[0];
    if (line[0]=='#' || line[0]==';' || line[0]=='!' || line[0]=='/')
      continue;
    while (isspace(*cp))
      *cp;
    if (*cp==0) continue;
    n++;
  }

  if (n>0) {
    for (np=1, p=0; ; np *= 2,p++) {
      if (np >= n) break;
    }
  }
  dprintf(0,"Found n=%d, p=%d\n",n,p);
  strclose(instr);
  return p;
}


int parse_team (
		string fname,        /* input     file name */
		int p,               /* input     depth of draw */
		int maxteam,         /* input     max # teams in list */
		Team *team,          /* input/o   pointer to list of teams (filled in upon output) */
		int maxplayer        /* input     max # players in list */
    )
{
    stream instr;
    char *cp, line[MAXLINELEN];
    string *words;
    int idx, len, nwords, rank=-1, count=0;
    int i, i1, i2, n = maxteam;

    if (fname==0 || *fname==0) {
      if (maxplayer==0) maxplayer=maxteam;
      warning("Inserting dummy %d players",maxplayer);
      for (i=0; i<maxplayer; i++) {
	sprintf(team[i].name,"Player-%d",i+1);
	team[count].id = i+1;
      }
      return maxplayer;
    }

    instr = stropen(fname,"r");
    while (fgets(line,MAXLINELEN,instr) != NULL) {
        len = strlen(line);
        /* patch line, skip comments */
        if (line[len-1]=='\n') line[len-1]='\0';
        if (line[0]=='#' || line[0]==';' || line[0]=='!' || line[0]=='/') {
            continue;
        }
        if (count >= maxteam) warning("Too many entries; maxteam=%d",maxteam);
        words = burststring(line," \t");
        nwords = xstrlen(words,sizeof(string))-1;
        if (nwords < 1) continue;

        idx = 0;
#if 0
        cp = words[idx];
        if (isdigit(*cp)) {
            team[count].id = atoi(cp);
            idx++;
        } else
            team[count].id = count+1;
#else
        team[count].id = count+1;
#endif

        cp = words[idx];
        if (isdigit(*cp)) {
            team[count].rank = atoi(cp);
            idx++;
        } else
            team[count].rank = rank--;

        cp = words[idx];
        if (isdigit(*cp)) {
            team[count].region = atoi(cp);
            idx++;
        } else
            team[count].region = 0;

        cp = words[idx];
        if (isdigit(*cp)) {
            team[count].slot = atoi(cp);
            idx++;
        } else
            team[count].slot = 0;               /*   <--- default !!! */

        strcpy(team[count].name,words[idx++]);
        for (;idx < nwords; idx++) {
            cp = words[idx];
            if (*cp=='#') break;
            strcat(team[count].name," ");
            strcat(team[count].name,cp);
        }
        dprintf(2,"Reading: %d %s\n",team[count].id, team[count].name);
	if (maxplayer > 0 && count == maxplayer) {
	  warning("Only using first %d players/teams",maxplayer);
	  break;
	}
        count++;
    }

    if (count>0)
      qsort(team,count,sizeof(Team),rank_cmp);
    
    return count;
        
}
#if 0
int rank_cmp(
    Team *i,
    Team *j
    )
{
    return j->rank - i->rank;
}

int slot_cmp(
    Team *i,
    Team *j
    )
{
    return i->slot - j->slot;
}
#else
int rank_cmp(
    const void *i,
    const void *j
    )
{
  Team *ii = (Team *) i;
  Team *jj = (Team *) j;
  return jj->rank - ii->rank;
}

int slot_cmp(
    const void *i,
    const void *j
    )
{
  Team *ii = (Team *) i;
  Team *jj = (Team *) j;
  return ii->slot - jj->slot;
}

#endif


void make_draw(
    int nteam,      /* number of teams (<= n) */
    Team *team,     /* point to all teams */
    int p,          /* level */
    int n,          /* 2**p, a bit redundant, but handy to keep around */
    int nrandom     /* if non-zero, randomize beyond this */
    )
{
  int start, todo, i, j, k;
  int filled[256];
  int pick[257], idx[257];
  bool Qibf = FALSE;
  
  for (i=0; i<nteam; i++)                 /* debugging: report */
    dprintf(1,"Read: %d %d %d %s\n",i+1,
	    team[i].id, team[i].rank, team[i].name);
  
  for (i=0; i<=n; i++) {                  /* tag array to keep track */
    filled[i] = 0;
    idx[i] = i;
  }
  
  for (i=0; i<n/2; i++) {                 /* set the slots that need to be filled */
    k = byes[p][n/2-i-1];                 /* count backwards to figure out how to draw the first half */
    pick[i] = k%2 ? k+1 : k-1;            /* even and odd are handled differently */
    pick[i+n/2] = byes[p][i];             /* bottom halves, which may be partially used */
  }
  
  dprintf(1,"Slot fill order: ");         /* debug output */
  for (i=0; i<n; i++)
    dprintf(1," %d",pick[i]);
  dprintf(1,"\n");

  /*  nrandom controls how the names from the list are indexed into their
      respective slots. E.g. the number 1 and 2 are normally not randomized,
      numbers 3 and 4 can be kept in that order if nrandom=5,
      numbers 5-8 also if nrandom=9 etc.etc.
      Normally IBF rules claim 3/4 are equal, so nrandom=3 is needed.
      However, 5-8 are randomized within that range, 9-16 also etc.etc.

      Shadow seeds need to be programmed properly still....
  */
  
  if (nrandom) {                          /* shuffle the index array in blocks */
    if (nrandom <= 3) {
      dprintf(0,"Randomizing 3/4\n");
      k = MIN(2,nteam-2);     shuffle(k,&idx[2]);
    }
    if (nrandom <= 5) {
      dprintf(0,"Randomizing 5/8\n");
      k = MIN(4,nteam-4);     shuffle(k,&idx[4]);
    }
    if (nrandom <= 9) {
      dprintf(0,"Randomizing 9/16\n");
      k = MIN(8,nteam-8);     shuffle(k,&idx[8]);
    }
    if (nrandom <= 17) {
      dprintf(0,"Randomizing 17/32\n");
      k = MIN(16,nteam-16);   shuffle(k,&idx[16]);
    }
    if (nrandom <= 33) {
      dprintf(0,"Randomizing 33/64\n");
      k = MIN(32,nteam-32);   shuffle(k,&idx[32]);
    }
  }

  dprintf(1,"Idx randomized order: ");         /* debug output */
  for (i=0; i<n; i++)
    dprintf(1," %d",idx[i]);
  dprintf(1,"\n");

  for (i=0; i<nteam; i++) {                         /* fill the slots */
    team[i].slot = pick[idx[i]];
  }
}

/*
 * here's the old terribly buggy version of make_draw
 */
void make_draw_old(
    int nteam,      /* number of teams (<= n) */
    Team *team,     /* point to all taems */
    int p,          /* level */
    int n,          /* 2**p, a bit redundant, but handy to keep around */
    int nrandom     /* if non-zero, randomize beyond this */
    )
{
    int start, todo, i, j, k, idx=0;
    int filled[256];
    int pick[257];
    bool Qibf = FALSE;

    for (i=0; i<nteam; i++)                 /* debugging: report */
        dprintf(1,"Read: %d %d %d %s\n",i+1,
                team[i].id, team[i].rank, team[i].name);

    for (i=0; i<=n; i++)                    /* tag array to keep track */
        filled[i] = 0;

    /* set the first 16 slots; note these are symmetric in both halves 
     * and  not like e.g. in IBF draws 
     * I'm now adding the IBF rules using a hardcoded Qibf = FALSE
     */
    pick[0]  = 1;                          
    pick[1]  = n;                          
    pick[2]  = (Qibf ?   n/4+1        :  n/2+1              );    /* 3/4 */
    pick[3]  = (Qibf ?   n+1-(n/4+1)  :  n/2                );
    pick[4]  = (Qibf ?   0            :  n/4+1              );    /* 5/8 */
    pick[5]  = (Qibf ?   0            :  n+1-(n/4+1)        );
    pick[6]  = (Qibf ?   0            :  n+1-n/4            );
    pick[7]  = (Qibf ?   0            :  n/4                );
    pick[8]  = (Qibf ?   0            :  n/8+1              );   /* 9/16 */
    pick[9]  = (Qibf ?   0            :  n+1-(n/8+1)        );
    pick[10] = (Qibf ?   0            :  n+1-(n/4+n/8)      );
    pick[11] = (Qibf ?   0            :  n/4+n/8            );
    pick[12] = (Qibf ?   0            :  n/4+n/8+1          );
    pick[13] = (Qibf ?   0            :  n+1-(n/4+n/8+1)    );
    pick[14] = (Qibf ?   0            :  n+1-n/8            );
    pick[15] = (Qibf ?   0            :  n/8                );

    if (nrandom) {
        if (nrandom <= 3) {
            dprintf(0,"Randomizing 3/4\n");
            shuffle(2,&pick[2]);
        }
        if (nrandom <= 5) {
            dprintf(0,"Randomizing 5/8\n");
            shuffle(4,&pick[4]);
        }
        if (nrandom <= 9) {
            dprintf(0,"Randomizing 9/16\n");
            shuffle(8,&pick[8]);
        }
    }


    while (idx < MIN(nteam,16)) {           /* fill in the first 16 */
        team[idx].slot = pick[idx];
        if (filled[team[idx].slot]) error("%d already filled",idx+1);
        filled[team[idx].slot] = 1;    
        idx++;
    }

    todo = n;
    for (i=0; i<n; i++)
    	if (filled[i+1]) todo--;
    dprintf(0,"After filling slots 1..16 there are %d left to do; nteam=%d\n",
		todo,nteam);
    
    while (idx < nteam) {    /*  Fill in final 16-8 "randomly" */
    	k = (int) xrandom(0.0, (double) todo) + 1;
    	dprintf(2,"Random: k=%d\n",k);


    	start=0;
	while (k>0) {
	    start++;
	    if (start>n) start=1;
	    if (!filled[start]) k--;
	}
	dprintf(2,"Hit at slot=%d filled=%d\n",start,filled[start]);

	team[idx].slot = start;
        if (filled[team[idx].slot]) error("Trailing already filled");
	filled[team[idx].slot] = 1;
	dprintf(1,"Putting %d into slot %d\n",idx+1,team[idx].slot);
	if (++idx >= nteam) return;
    	todo--;
    }
}

/*
 *  take an array idx[] of length n, and shuffle its elements
 *
 *  Note that n<=0 is allowed and nothing interesting will happen
 */

void shuffle(int n, int *idx)
{
    int k, l, tmp;

    for (k=n-1; k>0; k--) {
        l = (int)xrandom(0.0,(double)k+1);
        if (k==l) continue;
        tmp= idx[k];
        idx[k] = idx[l];
        idx[l] = tmp;
    }
}


void flip(int n, int *p)
{
    int i, tmp;

    for (i=0; i < n/2; i++) {
        tmp = p[n-1-i];
        p[n-1-i] = p[i];
        p[i] = tmp;
    }
            
}

/*
 * this elegant (?) algorithm codifies for arbitrary draw size
 * the draw locations where the byes are supposed to be located
 * according to IBF/USAB rules.
 * I still use the 1/n .... n/n notation, i.e. slots in the
 * draw are numbered 1..n, teams that populate these slots
 * are 1..m (actually called nteam in the code)
 *

 byes[0]  N/A
 byes[1]  2
 byes[2]  3 2
 byes[3]  5 4 7 2
 byes[4]  9 8 13 4 11 6 15 2
 byes[5]  17 16 25 8 21 12 29 4 19 14 27 6 23 10 31 2
 byes[6]  33 32 49 16 41 24 57 8 37 28 53 12 45 20 61 4 35 30 51 14 43 22 59 6 39 26 55 10 47 18 63 2

 */

void setup_byes(int pmax)
{
  int i, n, p;

  byes = (int **) allocate((pmax+1)*sizeof(int *));
  byes[0] = NULL;  // actually never need this, p > 0
  
  for (p=1, n=1; p<=pmax; p++, n*=2) {
    byes[p] = (int *) allocate((n+1)*sizeof(int));
    if (p==1) {          // first iteration, actually not useful
      byes[p][0] = 2;
      byes[p][1] = 0;    // terminator
    } else if (p==2) {   // second iteration, setup for the remainder
      byes[p][0] = 3;
      byes[p][1] = 2;
      byes[p][2] = 0;    // terminator
    } else {
      for (i=0; i<n/2; i++)     // recurse from previous p
	byes[p][i] = i%2 ? 2*byes[p-1][i] : 2*byes[p-1][i]-1;
      for (i=n/2; i<n; i++)     // recurve from top half of this p
      	byes[p][i] = i%2 ? byes[p][i-n/2]-2 : byes[p][i-n/2]+2;
      byes[p][n] = 0;
    }
    dprintf(1,"### p=%d n=%d:",p,n);
    for (i=0; byes[p][i] > 0; i++)
      dprintf(1," %d",byes[p][i]);
    dprintf(1,"\n");
  }
}

