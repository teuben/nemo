/*
 * TABLOVAS: Lovas spectral line catalogue browser
 *          See the man page for many more details.
 *
 *	29-apr-92  V0.0    toy version, cloned off tabpp                    PJT
 *			   ifdeffed as hash search table for the species  
 *			   which would speedup any wildcard mathing to be done
 * 
 */
/**************** NEMO * UNIX INCLUDE FILES **********************************/ 

#include <stdinc.h>
#include <getparam.h>
#include <math.h>
#include <hash.h>

/**************** COMMAND LINE PARAMETERS ************************************/

string defv[] = {
    "in=???\n       Input table",
    "species=all\n  Which species (note need for upper case)",
    "col=\n         List of columns to see [(name,freq)|all]",
    "freq=\n        Ranges in frequency to show, and units ([Mhz]|Ghz|THz)",
    "cutoff=\n      Cutoff in intensity *** ",
    "tab=\n         Output table [terminal](** not implemented **)",
    "select=\n      C-style selection criteria (** not implemented **)",
    "fmt=\n         Alternative format for column output (** not implemented **)",
    "mode=f\n       Extra output of all species?",
    "VERSION=1.1\n  8-may-92 PJT",
    NULL,
};

string usage = "Lovas spectral line catalogue browser";

/**************** GLOBAL VARIABLES ********************************************/

real freqmin=0.0,
     freqmax=999999999.999,
     factor=1.0,
     lfactor=0.0,
     cutoff=-999999.9;

int scount, lcount;


/****************************** START OF PROGRAM *****************************/

nemo_main()
{
    stream instr;
    string mode, tab, *burststring(), sconc();
    string *species, *cols, *freqs, names;
    real tmp;
    int  i;
    bool Qalls, Qallc, Qspecies, scanopt();
    char sname[128], sfile[128], *dname, *cp, *getenv();
    double log10();

    species = burststring(getparam("species"),", ");
    Qalls = (xstrlen(species,sizeof(string))<=1 || streq(species[0],"all"));

    cols = burststring(getparam("col"),", ");
    Qallc = (xstrlen(cols,sizeof(string))<=1 || streq(cols[0],"all"));

    if(hasvalue("cutoff")) cutoff = getdparam("cutoff");

    freqs = burststring(getparam("freq"),", ");
    switch(xstrlen(freqs,sizeof(string))-1) {
     default:
        warning("Bad parsing: assume whole freq range...\n");
     case 0:
        break;
     case 1:
        nemoinpd(freqs[0],&freqmin ,1);
        break;
     case 2:
        nemoinpd(freqs[0],&freqmin ,1);
        nemoinpd(freqs[1],&freqmax ,1);
        break;
     case 3:
        nemoinpd(freqs[0],&freqmin ,1);
        nemoinpd(freqs[1],&freqmax ,1);
        strlower(freqs[2]);
        if(streq(freqs[2],"mhz"))
            factor = 1.0;
        else if (streq(freqs[2],"ghz")) {
            warning("Units in GHz: note freq, etc. are modified");
            factor = 0.001;
        } else if (streq(freqs[2],"thz")) {
            warning("Units in GHz: note freq, etc. are modified");
	    factor = 0.000001;
        } else
            error("%s is an illegal frequency unit: must be (Mhz,Ghz,Thz)",
                  freqs[2]);
        break;
    }
    lfactor = log10(factor);
    if(freqmin>freqmax) {
          warning("Toggling freq's since they much be sorted");
          freqmin = tmp;
          freqmin = freqmax;
          freqmax = tmp;
    }
    Qspecies = getbparam("mode");


    check_line(getparam("in"), Qalls, Qallc, species, cols, Qspecies);
}

patchzero(name)
char *name;
{
    while (*name) if(*name == ' ') *name++ = '0'; else name++;
}

patchnull(name)
char *name;
{
    while (*name) if(*name == ' ') *name++ = '\0'; else name++;
}

/* MYSUBSTR: local substring copy */

char *mysubstr(line,n,l)
char *line;
int *n,l;
{
    local char word[100];

    strncpy(word,&line[*n],l);
    word[l] = NULL;
    *n += l;

    return(word);
}

/* STRLOWER: convert a string to lower case */

strlower(s)
char *s;
{
    while (*s) {
       if(isupper(*s)) *s = tolower(*s);
       s++;
    }
}
       

check_line(file, Qalls, Qallc, species, cols, Qspecies)
string file;
bool Qalls, Qallc, Qspecies;
string *species, *cols;
{
    stream instr = stropen(file,"r");
    char line[256], *mysubstr();
    char molecule[32], smol[32], *sm;
    bool found, wild();
    real freq;
    int  i, n, nspecies=0, ncols=0, nlines=0, nout=0, linelen, *p;
    struct Hash_Table *h;

    if (!Qalls) nspecies = xstrlen(species,sizeof(string))-1;
    if (!Qallc) ncols    = xstrlen(cols,sizeof(string))-1;

    h = init_Hash_Table();

    if (fgets(line,256,instr)==NULL)
        error("Reading first testline");
    switch(line[132]) {
      case '\n': dprintf(0,"LOVAS: Ascii version\n");
                 linelen=133; 
                 break;
      case ' ':  dprintf(0,"LOVAS: Binary version\n");
                 linelen=132; 
                 break;
      default:   linelen = strlen(line);
                 warning("Unexpected linelength %d used\n",linelen);
                 break;
    }
    rewind(instr);

    while( fgets(line,linelen+1,instr) ) {
	line[linelen] = '\0';
        if (strlen(line) != linelen) {
		warning("Skipping line[%d]: %s\n",strlen(line),line);
		continue;
	}
	nlines++;
        if(line[linelen-1]=='\n') line[linelen-1] = '\0';  /* patch ascii */
        
        n=1;                                        /* 1x */
        strcpy(molecule, mysubstr(line,&n,14));     /* A14 */
        freq = atof(mysubstr(line,&n,10));          /* F10.3 */

	/* change some units for output */
        freq *= factor;


        if(freqmin <= freq && freq <= freqmax) {
            if(nspecies>0) {
#if 0
                /* simple linear search */
                found = FALSE;
                for(i=0;i<nspecies;i++)
                    if(strncmp(molecule,species[i],strlen(species[i]))==0) {
                        found=TRUE;
                        break;
                    }
                if(!found) continue;    /* continue with next line */
#else
                /* search while hashing all species into lookup table */
                strcpy(smol,molecule);
                patchnull(smol);
                
                if( (p=(int *)get_hash(h,smol))==NULL) { /* new one: hash it */
                    p = (int *) allocate(sizeof(int));
                    *p = 0;
                    sm = (char *) allocate(strlen(smol)+1);
                    strcpy(sm,smol);
                    for(i=0; i<nspecies; i++)
                        if (wild(species[i],smol)) {
                            *p = i+1;
                            break;
                        }
                    if (!put_hash(h,sm,p)) 
			error("put_hash %s; idx=%d record=%d\n",sm,*p,nlines);
                    if (Qspecies) printf("Species: %s\n",sm);
                } 
                if(*p==0) continue;    
#endif
            } /* species-test */
            nout++;
            if(Qallc)
                printf("%s\n",line);
            else
                printf("%-14s %f\n",molecule,freq);
        } /* freq-test */

    } /* get-line */

    strclose(instr);
    dprintf(0,"Read %d lines; output %d lines\n",nlines, nout);
}

/*
 *  Wildcard match: currently only matches without '*' and '?' expansion
 *      The match is also considered false if the two strings are unequal
 *      in length, though the first blank is considered EOS also....
 */

#define EOS(a)   ( (a)==NULL || (a)==' ')

bool wild(a,b)     
char *a, *b;
{
    if (EOS(*a) || EOS(*b)) return FALSE;


    while (*a==*b) {
        if(EOS(*a)) return TRUE;
        a++;
        b++;
    }
    return FALSE;
}
