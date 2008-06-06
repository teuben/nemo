/*
 * BODYTRANS.C: functions to look up or compile a body transformation
 * and return a pointer to the function so defined.
 *
 *  =>  Experimental version using a slow interpreter (fie.c)
 *	when loadobj.c is not functional on your host
 *
 *	10-sep-90  Created, but never finished  -  Peter Teuben
 *	22-mar-91  worked a bit more on it	PJT
 *
 *	<<< STILL TESTING - FEEL FREE TO FINISH IT OFF>>>
 *	[got side tracked by loadobjDL.c which seems to work, sort of]
 */

#include <stdinc.h>
#include <getparam.h>
#include <filefn.h>
#include <loadobj.h>
#include <strlib.h>

/*
 * BTRTRANS, BTITRANS: map name or expression to real or integer
 * valued function, and return a pointer.
 */

rproc btrtrans(expr)
string expr;			/* name or C expression */
{
    proc bodytrans();

    return ((rproc) bodytrans("real", expr, NULL));
}

iproc btitrans(expr)
string expr;			/* name or C expression */
{
    proc bodytrans();

    return ((iproc) bodytrans("int", expr, NULL));
}



/*
 * BODYTRANS: local workhorse function, which (compiles and) loads
 * the specified transformation, and returns a pointer.
 */

local int  funcmpld = 0;	/* count of functions compiled */
local string *ebase=NULL;       /* point to array of allocated expr */


local proc bodytrans(type, expr, fname)
string type;			/* type of function to return */
string expr;			/* name or C expression */
string fname;			/* not used in this _fie version */
{
    proc result;
    char code[256], *cp;
    string *sp;
    int codelen, idx, n,l,slot;

    /* first the expression needs to be parsed and bodyvariables
       replaced by their index %n (n=1,2,...) */

    
    if ( (l = strlen(expr)) <= 0) return((proc )NULL);
    slot = -1;
    if (ebase) {
        for (sp=ebase; *sp; sp++) {
            if (strlen(*sp)==l) {
                if (streq(*sp,expr)) {
                    dprintf(0,">>> Found old expr %s\n",expr);
                    slot = (int) (sp-ebase);
                    break;
                }
            }
        }
    }
    if (slot < 0) {             /* no slot found - so allocate one */
        dprintf(0,">>> Entering new expr %s\n",expr);
        funcmpld++;
        if (ebase==NULL) 
            ebase = (string *) allocate(2*sizeof(string *));
        else
            ebase = (string *) reallocate(ebase,(funcmpld+1)*sizeof(string *));
        ebase[funcmpld-1] = scopy(expr);        /* save expression */
        ebase[funcmpld] = NULL;
        codelen = 0;        /* count length of code we have so far */
        cp = expr;          /* cp will point where in 'expr' we are parsing */
        while (*cp) {
            if (isbodyvar(cp,&idx,&n)) {
                sprintf(&code[codelen],"%%%d",idx+1);
                cp += n;
                codelen = strlen(code);
            } else {
                while (n--)
                    code[codelen++] = *cp++;
            }
        }
        code[codelen] = '\0';
        printf("Code = %s\n",code);
        inifie(code);
        savefie(funcmpld-1);
    } else
        loadfie(slot);

    stop(0);
}

local string bodyvar[] = {	/* Bodyvars must obey 'a' <= ch <= 'z' */
      "i","t","m",              /*   in order for the fast 'isbodyvar' to */
      "x","y","z",              /*   work */
      "vx","vy","vz",
      "phi",
      "ax","ay","az",
      "aux","key",
      NULL
};

/*
 *
 *  name (in)   start of string to check if a bodyvar
 *  idx  (out)  index in bodyvar array which bodyvar was found
 *  len  (out)  length of bodyvar string if one found
 *		 if none found length of name that can be skipped
 *		 because a bodyvar can only contains ['a'..'z']
 *
 *  isbodyvar (out)  :  0   is not a bodyvar
 *                      1   is a bodyvar
 */

local int isbodyvar(name,idx,len)
char *name;
int *idx, *len;
{
    int i, l;
    
    *len = 1;
    if (*name < 'a' || *name > 'z') return(0);

    for (i=0; bodyvar[i] != NULL; i++) {
        l = strlen(bodyvar[i]);
        if (strncmp(name,bodyvar[i],l)==0) {        /* got one ? */
            name++;                         /* check next char */
            if (*name >= 'a' && *name <= 'z') {   /* guess not */
                return(0);
            }
            *len = l;
            *idx = i;
            return(1);
        }
    }
    l = 0;  /* not found; see how much we can efficiently skip ... */
    while (*name && *name >= 'a' && *name <= 'z') {   /* skip alpha part */
        l++;
        name++;
    }
    while (*name && (*name < 'a' || *name > 'z')) {      /* skip non-alpha */
        l++;
        name++;
    }
    *len = l;
    return(0);
}



#include <bodytrans.h>

local float bin[15];   /* i,t,m,x,y,z,vx,vy,vz,phi,ax,ay,az,aux,key */
local float bout[1], berr=0.0;
local int   bn=1;	/* never change this ! */

local real _fnr_1(b,t,i)
Body *b;
real t;
int i;
{
    if (loadfie(1)<0) error("bodytrans_1\n");
    loadbpars(b,t,i);
    dofie(bin,&bn,bout,&berr);
    return ((real) bout[0]);
}

local loadbpars(b,t,i)      /* 3D only for now and standard Body */
Body *b;
real t;
int i;
{
    bin[0] = i;
    bin[1] = t;
    bin[2] = m;
    bin[3] = x;
    bin[4] = y;
    bin[5] = z;
    bin[6] = vx;
    bin[7] = vy;
    bin[8] = vz;
    bin[9] = phi;
    bin[10] = ax;
    bin[11] = ay;
    bin[12] = az;
    bin[13] = aux;
    bin[14] = key;
}


#ifdef TESTBED

/*
 * Test function, also useful when trying out new transformations.
 */
string defv[] = {
    "expr=x\n		Expression to test/save",
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
    "file=\n		Filename to save expression in",
    "btnames=\n		BTNAMES filename to regenerate .o files",
    "VERSION=1.3b\n	14-mar-90 PJT",
    NULL,
};

local getvparam();

main(argc, argv)
string argv[];
{
    string expr, type, fname;
    char *cp;
    Body b;
    real t;
    int i;
    rproc rtrans;
    iproc itrans;

    initparam(argv, defv);

    fname = getparam("btnames");
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
        printf("%s = %g \n", expr, (*rtrans)(&b,t,i));
    } else if (type[0] == 'i') {
        itrans = (iproc) bodytrans("int", expr, fname);
        printf("%s = %d \n", expr, (*itrans)(&b,t,i));
    } else
        dprintf(0,"Warning: not a valid type, must be real or int\n");
}

local getvparam(vec, name)
vector vec;
string name;
{
    string *burststring(), *cmpt;
    int i;
    real atof();

    cmpt = burststring(getparam(name), ", ");
    for (i = 0; i < NDIM; i++)
        vec[i] = atof(*cmpt != NULL ? *cmpt++ : "0.0");
}

#endif

