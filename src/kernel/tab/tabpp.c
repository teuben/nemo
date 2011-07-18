/*
 * TABPP: Poynter-Pickett spectral line catalogue browser
 *          See the man page for many more details.
 *
 *      29-feb-92  V1.0     A toy version created                  PJT
 *	 9-mar-92      a   bit more doc, allow Thz units
 *			   handle $ENVNAME in dir=		   PJT
 *	12-mar-92      c   added cutoff=			   PJT
 *      20-mar-92      d   warning when not using mhz - to appease joe
 *      24-mar-92  V1.1    added computation of \mu, the electric
 *                         dipole matrix element                   PJT
 *	 3-oct-92  V1.2    implemented select=                     PJT
 *	19-jan-95  V1.3    fixed up for the new (>1994) catalog    PJT
 *				(rely on zcat to do .Z and .gz)
 * 	11-sep-95  V1.3a   fixed parsing error $PP_DIR		   pjt
 *	 9-oct-95      b   new .cat extension, instead of .ed4	   pjt
 *      15-may-98  V1.4    merged two confusing version            pjt
 *      17-jul-11  V1.5    fixed copying into colval[]             pjt
 */
/**************** NEMO * UNIX INCLUDE FILES **********************************/ 

#include <stdinc.h>
#include <getparam.h>
#include <filefn.h>
#include <strlib.h>
#include <cgs.h>		/* CGS constants */

/**************** COMMAND LINE PARAMETERS *************************************/

string defv[] = {
    "species=\n                    List of species [all]",
    "mode=???\n                    Mode of output {species|lines}",
    "col=\n                        List of columns to see [all]\n\
        In species mode: (tag,name,nline,qlog,version)|all\n\
        In lines mode:   (freq,err,lgint,dr,elow,gup,tag,qnform,qn1,qn2)|all",
    "freq=\n    Ranges in frequency to show, and units ([Mhz]|Ghz|THz)",
    "cutoff=\n  Cutoff in LGINT (log of intensity)",
    "tab=\n     Output table [terminal](** not implemented **)",
    "select=\n  C-style selection criteria (** not implemented **)",
    "fmt=\n     Alternative format for column output (** not implemented **)",
    "dir=$PP_DIR,c%s.cat,catdir.cat\n  Dir, species-fmt, sp-dirfile",
    "VERSION=1.5\n               17-jul-11 PJT",
    NULL,
};

string usage = "Poynter-Pickett spectral line catalogue browser";

/**************** GLOBAL VARIABLES ********************************************/

real freqmin=0.0,
     freqmax=999999999.999,
     factor=1.0,
     lfactor=0.0,
     cutoff=-999999.9;

string fmt1 = "     FREQ         ERR     LGINT  DR     ELOW   GUP   TAG   NAME   MU^2\n";
string fmt2 = "%14.4f %8.4f %8.4f %2d %10.4f %3d %7d %s %f\n";
string fmt3 = "__________________________________________________________\n";
string fmt4 = "   NAME         TAG    Nlines\n";
string fmt5 = 
"   TAG NAME,         NLINE QLOG_1 QLOG_2 QLOG_3 QLOG_4 QLOG_5 QLOG_6 QLOG_7 VERSION\n";
string fmt6 =
"        FREQ     ERR    LGINT DR    ELOW   GUP   TAG QNFORM QN1(6)   QN2(6)\n";

string fmt7 = "%14.4f%8.4f%8.4f%s\n";

#define MAXSPECIES 350

struct species {
    char tag[7];
    char name[15];
    int nline;
    float qlog[7];
    int version;
} species[MAXSPECIES];

string colnam[] = {
    "freq", "err", "lgint", "dr", "elow", "gup", "tag", "qnform", "mu2",
    NULL};

string colval[9] ;
int scount, lcount;

string defdir = "/lma/spectra";
string deffmt = "c%s.ed4";
string defcat = "catdir.ed4";

extern string *burststring(string, string);
extern bool scanopt(string, string);

local string *burstcommand(string);
local char *mysubstr(string,int *, int);
local void strlower(char *s);

/****************************** START OF PROGRAM *****************************/

nemo_main()
{
    stream instr;
    string mode, tab;
    string *dirs, *usespecies, *cols, *freqs, *sels, names;
    real tmp;
    int  i;
    bool lmode, Qcolall, Qselall;
    char sname[128], sfile[128], *cp;
    string dname;

    tab = getparam("tab");
    mode = getparam("mode");
    lmode = (*mode=='l' || *mode=='L');

    names = getparam("species");
    usespecies = burststring(names,", ");
    cols = burststring(getparam("col"),", ");
    Qcolall = (xstrlen(cols,sizeof(string))<=1 || streq(cols[0],"all"));
    dirs = burststring(getparam("dir"),", ");
    sels = burststring(getparam("select"),",");   /* note only comma separated, not space */
    Qselall = (xstrlen(sels,sizeof(string))<=1);
    if(xstrlen(dirs,sizeof(string)) != 4)
        error("dir=$PP_DIR,sp-fmt,catname");
    if(hasvalue("cutoff")) cutoff = getdparam("cutoff");

    cp = dirs[0];
    if (*cp == '$') {
        cp++;
        if (*cp)
             dname = getenv(cp);
        else
             dname = getenv("PP_DIR");  /* If name given: try PP_DIR */
        if (dname==NULL) dname = sconc(".","/");    /* or else, use current */
    } else {
    	dname = dirs[0];
    }
    if (dname[strlen(dname)-1] != '/')	/* UNIX: make sure '/' is last char */
        dname = sconc(dname,"/");	/* for later catenation */

    freqs = burststring(getparam("freq"),", ");
    switch(xstrlen(freqs,sizeof(string))-1) {
     default:
        warning("Bad parsing: assume whole freq range...\n");
     case 0:
        break;
     case 1:
        nemoinpr(freqs[0],&freqmin ,1);
        break;
     case 2:
        nemoinpr(freqs[0],&freqmin ,1);
        nemoinpr(freqs[1],&freqmax ,1);
        break;
     case 3:
        nemoinpr(freqs[0],&freqmin ,1);
        nemoinpr(freqs[1],&freqmax ,1);
        strlower(freqs[2]);
        if(streq(freqs[2],"mhz"))
            factor = 1.0;
        else if (streq(freqs[2],"ghz")) {
            warning("Units in GHz: note freq, err and lgint are modified");
            factor = 0.001;
        } else if (streq(freqs[2],"thz")) {
            warning("Units in GHz: note freq, err and lgint are modified");
	    factor = 0.000001;
        } else
            error("%s is an illegal frequency unit: must be (Mhz,Ghz,Thz)",
                  freqs[2]);
        break;
    }
    lfactor = log10(factor);
    if(freqmin>freqmax) {
          warning("Toggling freq's since they must be sorted");
          freqmin = tmp;
          freqmin = freqmax;
          freqmax = tmp;
    }

    read_species(sconc(dname,dirs[2]), cols, lmode);


    if(lmode) {         /* mode=line ? */
        if (check_species(usespecies) == 0)
            error("No species left to search for");
        if(*names==NULL || streq(names,"all"))
            dprintf(1,"[Be patient, this operation may take a while]\n");
	if (Qcolall) {
            dprintf(0,fmt6);
	} else
	    dprintf(0,fmt1);
        for (i=0; i<scount; i++) {      /* loop over all species */
            sprintf(sname,dirs[1],species[i].tag);
            sprintf(sfile,"%s%s",dname,sname);
            if(*names && !streq(names,"all")) {
                if(!scanopt(names,species[i].name))
                    continue;
            }
            read_lines(sfile,&species[i],Qcolall,sels);
        }
    }
}

patchzero(char *name)
{
    while (*name) if(*name == ' ') *name++ = '0'; else name++;
}

patchnull(char *name)
{
    while (*name) if(*name == ' ') *name++ = '\0'; else name++;
}



/**********************************************************************/
read_species(fname, cols, lmode, Qall)		/* read data from input file */
string fname;
string *cols[];
bool lmode, Qall;
{
    char line[128], tag[6+1], name[14+1], note[1+1];
    int  nline, version, tag1, tag2;
    float qlog[7];
    stream instr;
    int i, n;

    instr = stropen(fname,"r");
    dprintf (2,"Reading speciesfile \n");
    if(!lmode) {
        if (Qall)
            dprintf(0,fmt5);
        else
            dprintf(0,fmt4);
    }

    scount = 0;
    while (fgets(line,128,instr) != NULL) {
        line[strlen(line)-1] = '\0';
        n = 0;
        strcpy(tag,mysubstr(line,&n,6));  tag[6] = '\0';
        patchzero(tag);  /* replace blanks by '0' */
        n++;
        strcpy(name,mysubstr(line,&n,14));  name[14] = '\0';
        patchnull(name);   /* replace blanks by '\0' */
        nline = atoi(mysubstr(line,&n,5));
        for(i=0; i<7; i++)
            qlog[i] = atof(mysubstr(line,&n,7));
        version = atoi(mysubstr(line,&n,2)); 
        strcpy(note,mysubstr(line,&n,1)); note[1] = '\0';

        strcpy(species[scount].tag,tag);
        strcpy(species[scount].name,name);
        species[scount].nline = nline;
        for(i=0; i<7; i++)
                    species[scount].qlog[i] = qlog[i];
        species[scount].version = version;

        scount++;
        lcount += nline;

        if (!lmode) {
            if (Qall)
                printf("%s\n",line);
            else
                printf("%-14s %-6s %5d\n",name,tag,nline);
        }
    }
    dprintf(0,"%s: %d species read, with a total of %d lines\n",
              fname,scount,lcount);
    strclose(instr);
}

read_lines(fname,spp,Qall,sels)	/* read species data from input file */
string fname;			/* pointed to by 'instr'     */
struct species *spp;            /* the species this applies to */
bool Qall;                      /* print all of one line? */
string *sels;                   /* list of selection criteria applied to each line */
{
    char line[128], zcatcmd[256], *cp;
    string ext;
    float freq, err, lgint, elow, eup, mu2;
    int   dr, gup, tag, qnform;
    int   qn1[6], qn2[6];
    int i, n;
    bool pipe;
    stream instr;

    dprintf (2,"Reading lines file %s\n",fname);
    ext = extension(fname);
    if (streq(ext,"Z") || streq(ext,"gz")) {
        pipe = TRUE;
        sprintf(zcatcmd,"zcat %s",fname);
        instr = popen(zcatcmd,"r");
    } else {
        pipe = FALSE;
        instr=stropen(fname,"r");
    }
    for (i=0; i<9; i++)  colval[i] = allocate(32);

    while (fgets(line,128,instr) != NULL) {          /* read row from line catalogue */
        line[strlen(line)-1] = '\0';
	dprintf(3,"%s\n",line);
        n = 0;
        /***  "freq", "err", "lgint", "dr", "elow", "gup", "tag", "qnform",   ***/
	/*    1041.4215   .0026 -8.7541 3   69.4458 23  46005 30311 3 8      11 3 9  */
        /*    1420.4058  0.0000 -9.0612 0    0.0000  3  -1001  22 1 1         1 0    */
	/*1234567890123........xxxxxxxx--xxxxxxxxxx---xxxxxxx----                    */

        cp = mysubstr(line,&n,13); freq =    atof(cp);  strcpy(colval[0],cp);
        cp = mysubstr(line,&n,8);  err  =    atof(cp);  strcpy(colval[1],cp);
        cp = mysubstr(line,&n,8);  lgint=    atof(cp);  strcpy(colval[2],cp);
        cp = mysubstr(line,&n,2);  dr   =    atoi(cp);  strcpy(colval[3],cp);
        cp = mysubstr(line,&n,10); elow =    atof(cp);  strcpy(colval[4],cp);
        cp = mysubstr(line,&n,3);  gup  =    atoi(cp);  strcpy(colval[5],cp);
        cp = mysubstr(line,&n,7);  tag  =    atoi(cp);  strcpy(colval[6],cp);
        cp = mysubstr(line,&n,4);  qnform  = atoi(cp);  strcpy(colval[7],cp);
        for(i=0; i<6; i++)
            qn1[i]  =  atoi(mysubstr(line,&n,2));
        for(i=0; i<6; i++)
            qn2[i]  =  atoi(mysubstr(line,&n,2));

	/* compute \mu when still in Mhz units */
        eup = elow + freq/29978.3;      /* E_upper in cm^-1 */
	mu2 = 1.502e11 * pow(10.0,(double)lgint+spp->qlog[0])
            * exp(eup/208.55) / (freq*freq*gup);
        sprintf(colval[8],"%g",mu2);

	/* change some units for output */
        freq *= factor;
        err  *= factor;
        lgint += lfactor;

        if(freqmin <= freq && freq <= freqmax && lgint>cutoff) {
            if(num_select(colnam, colval, sels))
                if(Qall)
                    printf(fmt7,freq,err,lgint,&line[29]);
                else
                    printf(fmt2,freq,err,lgint,dr,elow,gup,tag,spp->name,mu2);
        }
        dprintf(1,"%s",line);
    }
    if(pipe) pclose(instr);
    else   strclose(instr);
}

check_species(string *names)
{
    int i, n = xstrlen(names,sizeof(string))-1;
    bool found;
    string *sn;

    if (n==0) return(1);
    if (streq(names[0],"all")) return(1);

    for (sn=names; *sn; sn++) {
        found = FALSE;
        for(i=0; i<scount; i++) {
            if(streq(*sn,species[i].name)) {
                found = TRUE;
                break;
            }
        }
        if(!found) {
            warning("Species %s not present in catalogue\n",*sn);
            n--;
        }
    }
    return(n);
}



local char *mysubstr(char *line,int *n, int l)
{
    local char word[100];

    strncpy(word,&line[*n],l);
    word[l] = 0;
    *n += l;

    return word;
}

local string *burstcommand(string cmd)
{
    stream pstr;
    string *bc;
    char line[100];
    int n=0;

    bc = (string *) allocate((MAXSPECIES+1)*sizeof(string *));

    dprintf(3,"BURSTCOMMAND: CMD=%s\n",cmd);
    pstr = popen(cmd,"r");
    while (fgets(line,100,pstr) != NULL) {
        line[strlen(line)-1] = '\0';
        dprintf(3,"BURSTCOMMAND> %s\n",line);
        bc[n] = sconc(line,"");
        n++;
        if(n==MAXSPECIES) {
            warning("Too many species; MAXSPECIES=%d\n",MAXSPECIES);
            break;
        }
    }
    pclose(pstr);
    bc[n] = NULL;
    return(bc);
}

/* STRLOWER: convert a string to lower case */

void strlower(char *s)
{
    while (*s) {
       if(isupper(*s)) *s = tolower(*s);
       s++;
    }
}
       

