/*
 * BODYTRANS.C: functions to look up or compile a body transformation
 * and return a pointer to the function so defined.
 * 
 * public routines:
 *      rproc btrtrans(expr)
 *      iproc btitrans(expr)
 *
 *  -DTOOLBOX  version of this file can test bodytrans(5) files
 *  -DSAVE_OBJ will save bodytrans(5) files
 *  -DDLOPEN   assumes the dlopen() interface is used, sharable libraries
 *		need to be made instead. The .o will be a sharable
 *		object, or symbolically:
 *		(.c (cc[-pic) -> .o -> (mv) .z -> (ld]) -> .o )
 *
 *  History:
 *  xx-xxx-86  Created - Joshua Barnes
 *
 *  xx-mar-89  added a hash-expression parser when to save object files
 *          It works OK on BSD but clearly needs a better solution for 
 *          SYS5 where filenames are limited to 14 characters
 *
 *  xx-jul-89  bodytrans (a local function) has extra 'fname' parameter,
 *          to enhance power of TOOLBOX program
 *          Also split 'cc' into a 'cc' and 'mv', not using -o cc switch
 *          does not work on e.g. Ultrix (3.0) and SYS5
 *  
 *  19-sep-89  added ':' to hash-table (a bug really)
 *
 *  19-nov-89:  completed SYS5 port  (no leading underscores)
 *              and deleted hash technique, now using the local 
 *              XXX_bt functions, plus a simple file locking mechanism
 *              The -DSAVE_OBJ  option now uses the BTNAMES database.
 *
 *  14-mar-90:  V1.3 make GCC happy  [PJT]
 *
 *  10-sep-90:  V1.4 version for new DynamicLinker SUN OS 4.1 is promoting
 *			(dlopen() and friends)
 *		SUN manuals claim that when initialized data are used,
 *		not only .so is needed, also .sa - this is not implemented
 *		yet
 *  10-dec-90:  V1.4a - nemo_main() now
 *   5-nov-91   V1.4b - mapsys() was double defined in kernel/loadobj/mapsys.c 
 *                      already
 *   6-mar-92   Happy GCC2.0 - toolbox version uses hasvalue() now.
 *  13-mar-92   V1.5 include mathlinker to ensure loading math  functions
 *  15-nov-93   V1.6 extra check if requested filename length not too large
 *  10-jan-95   deleted atof() declaration for macro implementations
 *  13-dec-95   using new command ldso for dyn..obj.loader, 
 *              eliminated some nested extern's
 *  16-feb-97   attempt to support SINGLEPREC
 *   1-apr-01   NEMO V3 style .so file usage
 *
 *  Used environment variables (normally set through .cshrc/NEMORC files)
 *      NEMO        used in case NEMOOBJ was not available
 *      NEMOOBJ     normally points to $NEMO/obj/bodytrans
 *      BTRPATH     path of directories where to look for object files
 *	CFLAGS      if present, used in on-the-fly C compilation (only < V3)
 */

#include <stdinc.h>
#include <getparam.h>
#include <filefn.h>
#include <loadobj.h>
#include <ctype.h>
#include <mathlinker.h>		/* force inclusion of some basic math */

local proc   bodytrans();
local void   ini_bt(), end_bt(), make_bt();
local string get_bt(), put_bt();

/*
 * BTRTRANS, BTITRANS: map name or expression to real or integer
 * valued function, and return a pointer.
 */

rproc btrtrans(string expr)
{
    return (rproc) bodytrans("real", expr, NULL);
}

iproc btitrans(string expr)
{
    return (iproc) bodytrans("int", expr, NULL);
}

/*
 * DEFPATH: list of directories to search for body-trans functions,
 * used only if the environment variable BTRPATH is not set.
   !!!  this assumes there is /usr/nemo... and not $NEMO/... !!!
 */
 
#define DEFPATH  ".:/usr/nemo/obj/bodytrans"

/*
 * BODYTRANS: local workhorse function, which (compiles and) loads
 * the specified transformation, and returns a pointer.
 */

local bool havesyms = FALSE;       /* TRUE if symbols have been loaded */
local int  funcmpld = 0;                /* count of functions compiled */
#if defined(SAVE_OBJ)
local char edb[128], edbbak[128];     /* filenames of BTNAMES database */
local int  Qflock = 0;			/* used for file locking */
#endif

local proc bodytrans(type, expr, fname)
string type;                    /* type of function to return */
string expr;                    /* name or C expression */
string fname;                   /* optional filename for object file */
{
    char file[64], func[64], name[64], cmmd[256];
    char *sname, *cp, *cflags;
    string btrpath;
    string fullfile, hexpr;
    stream cdstr;
    proc result;

    if (! havesyms) {
        mysymbols(getargv0());
        ini_bt();
        havesyms = TRUE;
    }
    if (fname != NULL && *fname != 0)
        hexpr = fname;               /* set expression to 'filename' */
    else
        hexpr = expr;                  /* point to expr to try first */
    if ((int)strlen(hexpr) < 15) {
#if 1
        sprintf(file, "bt%c_%s.so", type[0], hexpr);     /* make filename */
#else
        sprintf(file, "bt%c_%s.o", type[0], hexpr);     /* make filename */
#endif
        dprintf(2,"bodytrans: trying file %s\n",file);
        btrpath = getenv("BTRPATH");
        fullfile = pathfind(btrpath != NULL ? btrpath : DEFPATH, file);
    } else {
        warning("bodytrans: skipping attempt to create file %s",hexpr);
        fullfile = NULL;
    }
        
    if (fullfile != NULL) {                  /* found a file already */
        loadobj(fullfile);                                /* load it */
        sprintf(func, "bt%c_%s", type[0], hexpr);    /* generic name */
        mapsys(func);                                     /* remap it */
    } else if ( (cp=get_bt(expr)) != NULL) {  /* use BTNAMES database */
	if (fname!=NULL && *fname!=0)
	    dprintf(0,"Found %s in %s\n",expr,cp);
	strcpy(file, getenv("NEMOOBJ"));             /* and found one */
	strcat(file, "/bodytrans/");
	strcat(file, cp);		 /* cp points to 'btX_NAME */
        strcat(file, ".o");
	loadobj(file);
        sprintf(func, "%s", cp);               /* generic symbol name */
        mapsys(func);                                     /* remap it */
    } else {                                           /* make a file */
        dprintf(0, "[bodytrans_new: invoking cc");
#if defined(SAVE_OBJ)
        dprintf(0, " +saving .o]\n");
#else
        dprintf(0, "]\n");
#endif
        if (fname!=NULL && *fname!=0)
            sprintf(name,"bt%c_%s",type[0],fname);      /* base name */
        else
            sprintf(name, "bt%c_%d", type[0], getpid()); /* temp name */
        dprintf(2,"bodytrans: base name = %s\n",name);

	sname = put_bt(expr, type[0], fname);	/* saved function name */
	sprintf(func,  "%s", sname);            /* general symbol name */
        mapsys(func);                                      /* remap it */
        dprintf(2,"bodytrans: func = %s\n",func);
        sprintf(file, "/tmp/%s.c", name);
        cdstr = fopen(file, "w");
        fprintf(cdstr, "#include <bodytrans.h>\n\n");
        fprintf(cdstr, "%s %s(b,t,i)\n", type, sname);/* use generic name */
        fprintf(cdstr, "Body *b;\nreal t;\nint i;\n");
        fprintf(cdstr, "{\n    return (%s);\n}\n", expr);
        fclose(cdstr);
	cflags = getenv("CFLAGS");
#if 1
        sprintf(cmmd, "cd /tmp;make -f $NEMOLIB/Makefile.lib %s.so",name);
#else
        sprintf(cmmd, "cd /tmp;cc %s -c %s.c",
		(cflags==NULL) ? "" : cflags,name);
#endif
        dprintf(2,"bodytrans: %s\n",cmmd);
        if (system(cmmd) != 0) {
#if defined(SAVE_OBJ)
            sprintf(cmmd,"rm -f %s",edbbak);    /* end file locking */
            system(cmmd);
#endif
            error("bodytrans(): could not compile\n");
	}
#if 1
        sprintf(file, "/tmp/%s.so", name);
#else
	sprintf(cmmd,"ldso /tmp/%s",name);
	dprintf(2,"bodytrans: %s\n",cmmd);
	if (system(cmmd) != 0)
		error("bodytrans: could not move link files\n");
        sprintf(file, "/tmp/%s.o", name);
#endif
        loadobj(file);

#if defined(SAVE_OBJ)
        if (!Qflock) {      /* copy when no file locking encountered */
#if 0
            sprintf(cmmd, 
            "cp /tmp/%s.so $NEMOOBJ/bodytrans/%s.so;chmod a+w %s",
             name,sname,edbbak);
#else
            sprintf(cmmd, 
            "cp /tmp/%s.o $NEMOOBJ/bodytrans/%s.o;chmod a+w %s",
             name,sname,edbbak);
#endif
            if (system(cmmd) != 0) {
                sprintf(cmmd,"rm -f %s",edbbak);    /* end file locking */
                system(cmmd);
                error("bodytrans(): could not copy\n");
            } else
                end_bt();       /* end the file locking */
        }
#endif            
    }
    result = findfn(func);
    if (result == NULL)
        error("Cant find %s (from findfn)", func);
    return result;
}

/*  
 * INI_BT: initialize some filenames for subsequent _BT functions
 *
 */

local void ini_bt()
{
#if defined(SAVE_OBJ)
    char *cp;

    cp = getenv("NEMOOBJ");
    if (cp!=NULL) {
        strcpy(edb, cp);
    } else {
        cp = getenv("NEMO");
        strcpy(edb, cp);
        strcat(edb, "/obj");
    } 
    strcat(edb, "/bodytrans/BTNAMES");
    strcpy(edbbak,edb);
    strcat(edbbak,".bak");
#endif
}

/*
 *  GET_BT:  find if there is a btX_NAME function in the BTNAMES 
 *           database with EXPR as expression. This is mainly useful 
 *           for expressions which do not translate into a proper
 *           filename, e.g. x+y
 */

#define NL 128

local string get_bt(expr)
string expr;
{
#if defined(SAVE_OBJ)
    int  nexpr;
    char *cp;
    static char line[NL];   /* warning: space may be used by pointer */
    FILE *fp;

    fp = fopen(edb, "r");
    if (fp==NULL)
        return(NULL);
    nexpr = strlen(expr);
    while (fgets(line, NL, fp))
        if (strncmp(expr,line,nexpr)==0) {
            cp = &line[nexpr];
            if (*cp != ' ')
                continue;                   /* not found yet */
            while (*cp==' ' || *cp=='\t')
                cp++;
            cp[strlen(cp)-1] = '\0';        /* strip '\n' */
            return(cp);
        }
#endif
    return(NULL);    			/* return as if nothing was found */
}

/* 
 * PUT_BT:  put away a compiled btX_NAME function in the database
 *	input:	expr	the full expresssion, e.g. x+y+2*x
 *              type    'r' for real, 'i' for int expression
 *		file	name to use for file, can be either a name
 *                      as in 'btX_NAME', or in case file==NULL
 *                      the type ('r' or 'i') is used to generate a
 *                      name and add it to the BTNAMES database
 *              A simple file locking mechanism has been used by
 *              creating a backup file first - the presence of a
 *              backup file signals file-locking of the master copy.
 */
 
#define MAXTRY 6

local string put_bt(expr,type,file)
string expr, file;
char type;
{
    FILE *fpi, *fpo;
    int  i;
    char *cp;
    static char line[NL];   /* warning: this space used by pointers */

#if defined(SAVE_OBJ)
    for (i=0; i<MAXTRY ; i++) {
        fpo = fopen(edbbak, "r");
        if (fpo==NULL) {
            fpo = fopen(edbbak, "w");
            if (fpo==NULL)
                error("Cannot open %s\n",edbbak);
            else
                break;
        } else {
            dprintf(1,"Waiting for someone to remove %s\n",edbbak);
            sleep(1);
        }
    }
    if (fpo!=NULL && i==MAXTRY) {
        sprintf(line,"ls -l %s",edbbak);
        system(line);
	dprintf(0,"Warning: bodytrans(5) expression %s not saved\n",expr);
	Qflock = 1;	/* flag that file was still locked by someone */
        sprintf(line,"bt%c__%d",type,++funcmpld);
        return(line);
    } else
        Qflock = 0;
        
    fpi = fopen(edb,"r");
    if (fpi!=NULL) {
	i=0;
        while(fgets(line,NL,fpi)) {
            fputs(line,fpo);
            i++;
        }
        fclose(fpi);
    } else
        dprintf(0,"[Creating %s\n",edb);
    if (file!=NULL && *file != 0)
        sprintf(line,"%s bt%c_%s\n",expr,type, file);
    else
        sprintf(line,"%s bt%c__%d\n",expr,type,i);
    fputs(line,fpo);
    fclose(fpo);

    line[strlen(line)-1] = '\0';    /* get rid of newline */
    cp = &line[strlen(expr)+1];     /* point to 'btX_NAME */
    return(cp);
#else
    sprintf(line,"bt%c__%d",type,++funcmpld);
    return(line);
#endif
}

/*
 * END_BT:  rename the backup file back to the original,
 *          which also signals the end of our simulated
 *          file-locking:
 *          as long as the backup exists, no one can work 
 *          on the BTNAMES database
 */

local void end_bt()
{
#if defined(SAVE_OBJ)
    if (!Qflock) {
       unlink(edb);
       link(edbbak,edb);
       unlink(edbbak);
    }
#endif
}

/*
 *  MAKE_BT:    (re)create all btX_NAME object files from this ascii 
 *              list in btname    *** NEVER USED ***
 */

local void make_bt(btname)
char *btname;
{
#if defined(SAVE_OBJ)
    FILE *fp;
    char *ep, *np, *tp, line[NL];
    int  i=0;

    fp = fopen(btname,"r");
    if (fp==NULL) {
        dprintf(0,"Could not open file %s\n",btname);
        return;
    }
    while (fgets(line,NL,fp)) {
        dprintf(1,"MAKE_BT: %s",line);
        ep= &line[0];
        np = ep;
        while (!isspace(*np))       /* skip the expression */
            np++;
        *np = '\0';                 /* terminate for 'ep' */
        np++;                       /* get to next char */
        while (isspace(*np))        /* skip spaces */
            np++;
        np++; np++;                 /* skip 2 ('b' and 't') */
        tp = np;                    /* 'i' or 'r' */
        np++; np++;                 /* skip 2 (type and '_') */
        np[strlen(np)-1] = '\0';    /* get rid of trailing '\n' */
        if (*np=='_')               /* renumber */
            sprintf(np+1,"%d",i);
        if (*tp == 'r')
            (void) bodytrans("real", ep, np);
        else if (*tp == 'i')
            (void) bodytrans("int", ep, np);
        else
            dprintf(1,"bodytrans: ignoring BTNAMES line (%s)\n",line);
        i++;
    }
#else
    dprintf(0,"bodytrans: not compiled with -DSAVE_OBJ\n");
#endif
}       

#ifdef TOOLBOX

/*
 * Test function, also useful when trying out new transformations.
 */

#include <vectmath.h>
#include <snapshot/body.h>

string defv[] = {
    "expr=x\n		Expression to test and/or save",
    "type=real\n	Function type : real/integer",
    "mass=0.125\n	A mass",
    "pos=0.5,-0.5,1.0\n	Positions",
    "vel=0.1,0.2,-0.3\n	Velocities",
    "phi=-1.2\n		Potential",
    "acc=-0.5,0.7,0.1\n	Accellerations",
    "aux=3.141592\n	Auxiliary",
    "key=3\n		Key",
    "t=2.5\n		Time",
    "i=1\n		Index",
    "file=\n		Filename to save expression in (btTYPE_FILE)",
    "btnames=\n		BTNAMES filename to regenerate .o files",
    "VERSION=2.0\n	1-apr-01 PJT",
    NULL,
};

string usage = "bodytrans manipulator";

local void getvparam();

void
nemo_main()
{
    string expr, type, fname;
    char *cp;
    Body b;
    real t;
    int i;
    rproc rtrans;
    iproc itrans;

    if (hasvalue("btnames")) {
        make_bt(getparam("btnames"));
        stop(0);
    }
    expr = getparam("expr");
    type = getparam("type");
    Mass(&b) = getdparam("mass");
    getvparam(Pos(&b), "pos");
    getvparam(Vel(&b), "vel");
    Phi(&b) = getdparam("phi");
    getvparam(Acc(&b), "acc");
    Aux(&b) = getdparam("aux");
    Key(&b) = getiparam("key");
    t = getdparam("t");
    i = getiparam("i");
    fname = getparam("file");

    if (type[0] == 'r') {
        rtrans = (rproc) bodytrans("real", expr, fname);
        cp = get_bt(expr);
        printf("%s = %g (%s)\n", expr, (*rtrans)(&b,t,i), (cp)?(cp):(""));
    } else if (type[0] == 'i') {
        itrans = (iproc) bodytrans("int", expr, fname);
        cp = get_bt(expr);
        printf("%s = %d (%s)\n", expr, (*itrans)(&b,t,i), (cp)?(cp):(""));
    } else
        dprintf(0,"Warning: not a valid type, must be real or int\n");
}

local void getvparam(vec, name)
vector vec;
string name;
{
    string *burststring(), *cmpt;
    int i;

    cmpt = burststring(getparam(name), ", ");
    for (i = 0; i < NDIM; i++)
        vec[i] = atof(*cmpt != NULL ? *cmpt++ : "0.0");
}

#endif
