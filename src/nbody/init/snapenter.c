/* snapenter.c - main, writesnap, getstring, ordinal */

/*
 *  snapenter.c:  for interactively entering an nbody snapshot
 *
 *    June 1988  -  Piet Hut  \@ Inst. f. Adv. Study, Princeton, NJ 08540, USA
 *
 *	10-jul-89  V2.0   updated for new filestruct()		PJT
 *			  does not opt to correct to Center of Mass coord's
 *	21-nov-90  V2.1   void/int repair - changed defaults/names  PJT
 *                        and fixed order of things done....
 *	24-mar-94  V2.2   added verbose=
 *       1-apr-97     a   fixed r-format bug
 *       1-apr-01     b   compiler warnings
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <extstring.h>

#include <filestruct.h>
#include <snapshot/snapshot.h>
#include <snapshot/body.h>
#include <snapshot/put_snap.c>

#include <ctype.h>

string  defv[] = {              /* DEFAULT INPUT PARAMETERS */
    "out=\n         Output file name",
    "nbody=0\n      Number of particles",
    "time=0.0\n     Time at which snapshot taken",
    "headline=\n    Verbiage for output",
    "verbose=t\n    Be verbose when input requested",
    "VERSION=2.2c\n 8-sep-01 PJT",
    NULL,
};

string usage = "enter an nbody snapshot interactively";

local string headline;                /* random text message */
local int    nobj;                    /* number of bodies in snapshot */
local real   tsnap;                   /* associated time of snapshot */
local Body   *btab = NULL;            /* pointer to Nbody snapshot */

local writesnap(stream outpt);
local string  getstring(void);
local string  ordinal(int i);

#if defined(SINGLEPREC)
local string rfmt = "%f";
#else
local string rfmt = "%lf";
#endif

/*-----------------------------------------------------------------------------
 *  nemo_main  --  interactively enters an nbody snapshot
 *            prompts for:  output file name (if not provided in argv)
 *                          number of particles (if not provided in argv)
 *                          particle data (mass, position and velocity
 *                            components for each particle).
 *            bugs:  only for 3 dimensions and Carthesian coordinates
 *-----------------------------------------------------------------------------
 */

void nemo_main()
{
    char    c;
    string oname;
    bool    Qverbose;
    stream  outstr;
    int     i;
    Body    *bp;

    oname = getparam("out");
    nobj = getiparam("nbody");
    tsnap = getdparam("time");
    headline = getparam("headline");
    Qverbose = getbparam("verbose");

    if (*oname == 0) {
        if (Qverbose) printf(">> output file name = ");
        oname = getstring();
    }
    if (strlen(oname) == 0)
        error("no output file name provided");
    outstr = stropen(oname,"w");

    if (nobj == 0) {
        if (Qverbose) printf(">> N = ");
        scanf("%d", &nobj);
    } 
    if (nobj < 1)
        error("nbody = %d , while nobj > 0 is required", nobj);

    btab = (Body *) allocate (sizeof(Body) * nobj);

    for (i=0, bp=btab; i<nobj; i++, bp++) {
        if (Qverbose) printf("## enter the data for the %d%s particle:\n",
                                                 i+1, ordinal(i+1));
        if (Qverbose) printf(">> mass = ");
        scanf(rfmt, &Mass(bp));
        if (Qverbose) printf(">> position:  x component = ");
        scanf(rfmt, &Pos(bp)[0]);
        if (Qverbose) printf(">>            y component = ");
        scanf(rfmt, &Pos(bp)[1]);  
        if (Qverbose) printf(">>            z component = ");
        scanf(rfmt, &Pos(bp)[2]);  
        if (Qverbose) printf(">> velocity:  x component = ");
        scanf(rfmt, &Vel(bp)[0]);  
        if (Qverbose) printf(">>            y component = ");
        scanf(rfmt, &Vel(bp)[1]);  
        if (Qverbose) printf(">>            z component = ");
        scanf(rfmt, &Vel(bp)[2]);  
        }
    warning("Center of mass not necessarily at the origin");
    writesnap(outstr);
    strclose(outstr);
    printf("Snapshot with %d particle(s) written to file %s\n",nobj,oname);
}

/*-----------------------------------------------------------------------------
 *  writesnap  --  writes a snapshot to an output file
 *                 accepts:  name, the name of the output file
 *-----------------------------------------------------------------------------
 */
local writesnap(stream outpt)
{
    int     bits = TimeBit | MassBit | PhaseSpaceBit;

    put_history(outpt);
    if (*headline != 0)                      /* non-trivial headline? */
        put_string(outpt, HeadlineTag, headline);
    put_snap(outpt,&btab,&nobj,&tsnap,&bits);
}

/*-----------------------------------------------------------------------------
 *  getstring  --  reads in a string from standard input,
 *                 terminated at the appearance of the first blank character
 *-----------------------------------------------------------------------------
 */
#define  BUFFERLENGTH   64

local string  getstring(void)
{
    char  buf[BUFFERLENGTH], *bp, c;
    int  i = BUFFERLENGTH;
    
    bp = buf;
    while (--i && !isspace(c = getchar()))
        *bp++ = c;
    *bp = 0;  /* terminate string */
    if (i)
        return (copxstr(buf, sizeof(char)));
    else
        error("getstring: string length exceeds %d characters\n",
                                                BUFFERLENGTH - 1);
}

/*-----------------------------------------------------------------------------
 *  ordinal  --  supplies the suffix to convert a cardinal number into
 *               an ordinal number
 *               accepts:  i, an integer
 *               returns:  a two-character string containing the suffix
 *               examples:  1  ==>   1st  
 *                         12  ==>  12th
 *                        123  ==> 123rd
 *-----------------------------------------------------------------------------
 */
local string  ordinal(int i)
{
    if ((i/10)%10 == 1)   /* all numbers ending with the last two digits */
        return("th");     /*   between 10 and 20 receive "th";           */
    else                  /* for all other numbers only the last digit   */
        switch (i%10)     /*   determines the suffix:                    */
            {
            case 1:  return("st"); 
                     break;
            case 2:  return("nd");
                     break;
            case 3:  return("rd");
                     break;
            default: return("th");
                     break;
            }
}

/* endof: snapenter.c */
