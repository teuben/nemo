/****************************************************************************/
/* GETPARAM.C: command-line processing functions.                           */
/* Copyright (c) 2001 by Joshua E. Barnes, Honolulu, Hawai`i.               */
/****************************************************************************/

#include "stdinc.h"
#include "getparam.h"
#include <string.h>

/*
 * PARAM: structure encoding parameter name and value.
 */

typedef struct {
    string name;                                /* name of parameter        */
    string value;                               /* value of parameter       */
    string comment;                             /* documentation string     */
    int flags;                                  /* for various options      */
} param;

/*
 * Local routines and definitions.
 */

local int countdefaults(string *);              /* number of defaults       */
local void setprogram(param *, string);         /* set 0th parameter        */
local void copydefaults(param *, string *);     /* set default parameters   */
local void checkhelp(param *, string);          /* help processing          */
local void printitem(string, string);           /* print param and comment  */
local void setarguments(param *, string *);     /* set command parameters   */
local void reqarguments(param *);               /* check manditory args     */
local param *findparam(string, param *);        /* look up parameter value  */
local string parname(string);                   /* extract param name       */
local string parvalue(string);                  /* extract param value      */

local param *paramvec = NULL;                   /* vector of parameters     */

local string progname = NULL;                   /* program name, for errors */

/*
 * INITPARAM: initalize parameter lists and handle special requests.
 */

void initparam(string *argv, string *defv)
{
    int nparam;
    param *pvec;

    progname = argv[0];                         /* initialize program name  */
    nparam = 1 + countdefaults(defv);           /* include argv0 in count   */
    pvec = (param *) allocate(sizeof(param) * (nparam + 1));
    setprogram(pvec, argv[0]);                  /* install 0th argument     */
    copydefaults(pvec, defv);                   /* set up default values    */
    checkhelp(pvec, argv[1]);                   /* give help if requested   */
    setarguments(pvec, argv);                   /* args override defaults   */
    reqarguments(pvec);                         /* complain if args missing */
    paramvec = pvec;                            /* install parameter vector */
}

/*
 * COUNTDEFAULTS: count number of default parameters.
 */

local int countdefaults(string *defv)
{
    int ndefault;
    string *dp;

    ndefault = 0;
    for (dp = defv; *dp != NULL; dp++)          /* loop over all defaults   */
        if (**dp != ';')                        /* if not a comment         */
            ndefault++;                         /* then count one more      */
    return (ndefault);
}

/*
 * SETPROGRAM: initialize the program name as parameter "argv0".
 */

local void setprogram(param *pvec, string argv0)
{
    pvec->name = "argv0";                       /* install 0th parameter    */
    pvec->value = argv0;                        /* set name from argv[0]    */
    pvec->comment = NULL;                       /* no comment for now       */
    pvec->flags = ARGPARAM;                     /* so user can't reset it   */
}

/*
 * COPYDEFAULTS: install default parameters and comments.
 */

local void copydefaults(param *pvec, string *defv)
{
    param *pp;
    string *dp, name, value;

    pp = pvec;                                  /* start with 0th param     */
    for (dp = defv; *dp != NULL; dp++)          /* loop over the defaults   */
        if (**dp == ';') {                      /* is this is a comment?    */
            if (pp->comment != NULL)            /* already have a comment?  */
                error("copydefaults: cant append comments\n");
            pp->comment = strdup(*dp + 1);      /* store comment in field   */
        } else {                                /* else its not a comment   */
            pp++;                               /* so move onto new param   */
            name = parname(*dp);                /* extract parameter name   */
            value = parvalue(*dp);              /* and parameter value      */
            if (name == NULL || value == NULL)  /* is either one missing?   */
                error("copydefaults: bad parameter %s\n", *dp);
            pp->name = strdup(name);            /* assign parameter name    */
            pp->value = strdup(value);          /* and parameter value      */
            pp->comment = NULL;                 /* clear comment field      */
            pp->flags = DEFPARAM;               /* set source to default    */
            if (streq(pp->value, "???"))        /* is this param required?  */
                pp->flags |= REQPARAM;          /* set the required flag    */
            if (**dp == '<')                    /* an input parameter?      */
                pp->flags |= INPARAM;           /* set the input flag       */
            else if (**dp == '>')               /* or an output parameter?  */
                pp->flags |= OUTPARAM;          /* set the output flag      */
            if (name[0] == '.')                 /* a hidden parameter?      */
                pp->flags |= HIDPARAM;          /* set the hidden flag      */
        }
    pp++;                                       /* past last real param     */
    pp->name = NULL;                            /* end list of parameters   */
}

/*
 * CHECKHELP: if requested, print out help mesaages and exit.
 */

local void checkhelp(param *pvec, string argv1)
{
    param *pp;
    char buf[128];

    if (argv1 != NULL && streq(argv1, "-clue")) {
                                                /* print brief help message */
        printf("%s", pvec->value);
        for (pp = pvec+1; pp->name != NULL; pp++)
            printf(" %s=%s", pp->name, pp->value);
        printf("\n");
        exit(0);
    }
    if (argv1 != NULL && streq(argv1, "-help")) {
                                                /* print full help message  */
       printitem(pvec->value, pvec->comment);
        for (pp = pvec+1; pp->name != NULL; pp++) {
            sprintf(buf, "  %s=%s", pp->name, pp->value);
            printitem(buf, pp->comment);
        }
        exit(0);
    }
}

local void printitem(string item, string comment)
{
    if (comment == NULL)
        printf("%s\n", item);
    else
        if (strlen(item) < 32)
            printf("%-32s  %s\n", item, comment);
        else
            printf("%s\n\t\t\t\t  %s\n", item, comment);
}

/*
 * SETARGUMENTS: override defaults with commandline arguments.
 */

local void setarguments(param *pvec, string *argv)
{
    bool scanpos;
    param *pp;
    string *ap, name;

    scanpos = TRUE;                             /* start scan by position   */
    pp = pvec;                                  /* start at 0th param       */
    for (ap = argv + 1; *ap != NULL; ap++) {    /* loop over command args   */
        name = parname(*ap);                    /* get param name, if any   */
        scanpos = scanpos && (name == NULL);    /* see how to match args    */
        if (scanpos) {                          /* matching by position?    */
            pp++;                               /* then move to next param  */
            if (pp->name == NULL)               /* make sure it exists      */
                error("%s: too many arguments\n", progname);
            pp->value = strdup(*ap);            /* OK, set new param value  */
        } else {                                /* else matching by name?   */
            if (name == NULL)                   /* make sure name was given */
                error("%s: nameless arg %s\n", progname, *ap);
            pp = findparam(name, pvec);         /* look for named param     */
            if (pp == NULL)                     /* make sure param exists   */
                error("%s: parameter %s unknown\n", progname, name);
            if (pp->flags & ARGPARAM)           /* must not already be set  */
                error("%s: parameter %s duplicated\n", progname, name);
            pp->value = strdup(parvalue(*ap));  /* OK, set new param value  */
        }
        pp->flags = (pp->flags & ~DEFPARAM) | ARGPARAM;
                                                /* switch source flag       */
    }
}

/*
 * REQARGUMENTS: print out short message on use.
 */

local void reqarguments(param *pvec)
{
    bool needarg;
    param *pp;

    needarg = FALSE;                            /* see if any args left out */
    for (pp = pvec+1; pp->name != NULL; pp++)   /* scan list of parameters  */
        if ((pp->flags & REQPARAM) &&           /* a required parameter?    */
              (pp->flags & DEFPARAM))           /* and still defaulted?     */
            needarg = TRUE;                     /* note missing args exist  */
    if (needarg) {                              /* list required arguments  */
        eprintf("Usage: %s", progname);
        for (pp = pvec+1; pp->name != NULL; pp++)
            if ((pp->flags & REQPARAM))
                eprintf(" %s=???", pp->name);
        error("%s: required arguments missing\n", progname);
    }
}

/*
 * GETPARAM: return value of parameter.  As a special case, requests for
 * "argv0" can be handled before paramvec has been set, to allow error
 * handling during the initialization process.
 */

string getparam(string name)
{
    param *par;

    if (paramvec == NULL)
        if (streq(name, "argv0") && progname != NULL)
            return (progname);
        else
            error("getparam: called before initparam\n");
    par = findparam(name, paramvec);
    if (par == NULL)
        error("getparam in %s: parameter %s unknown\n", progname, name);
    return (par->value);
}

/*
 * GETPARAMSTAT: return parameter flags, or zero if no such parameter
 * (note that all defined parameters have at least one bit set).
 */

int getparamstat(string name)
{
    param *par;

    par = findparam(name, paramvec);
    return (par != NULL ? par->flags : 0);
}

/*
 * GETIPARAM, GETDPARAM, GETBPARAM: get int, double, or bool parameter.
 */

int getiparam(string name)
{
    int val;
    string end;

    val = strtol(getparam(name), &end, 0);
    if (*end == 'k' || *end == 'K')
        return (1024 * val);
    if (*end == 'm' || *end == 'M')
        return (1024 * 1024 * val);
    return (val);
}

double getdparam(string name)
{
    return (atof(getparam(name)));              /* convert value to double  */
}

bool getbparam(string name)
{
    char *val;

    val = getparam(name);                       /* obtain value of param    */
    if (strchr("tTyY1", *val) != NULL)          /* is value true?           */
        return (TRUE);
    if (strchr("fFnN0", *val) != NULL)          /* is value false?          */
        return (FALSE);
    error("getbparam: %s=%s not bool\n", name, val);
    return (FALSE);                             /* keep compiler happy...   */
}

/*
 * FINDPARAM: look for named parameter in list of parameters.
 */

local param *findparam(string name, param *pvec)
{
    param *pp;

    for (pp = pvec; pp->name != NULL; pp++)
        if (streq(name, pp->name))
            return (pp);
    return (NULL);
}

/*
 * PARNAME: extract name from name=value string.
 * W A R N I N G :  returns ptr to static storage.
 */

local string parname(string arg)
{
    char *ap, *ep;
    static char namebuf[64];

    ap = (char *) arg;
    if (*ap == '<' || *ap == '>')
        ap++;
    strncpy(namebuf, ap, 63);
    namebuf[63] = NULL;
    ep = strchr(namebuf, '=');
    if (ep == NULL)                             /* not of form name=value?  */
        return (NULL);
    *ep = NULL;
    return (namebuf);
}

/*
 * PARVALUE: extract value from name=value string.
 */

local string parvalue(string arg)
{
    char *ep;

    ep = strchr(arg, '=');
    if (ep == NULL)
        return (NULL);
    return (ep + 1);
}
