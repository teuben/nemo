/* 
 *  FTOC:
 *      1)  extract a defv[] string array for a nemo_main()
 *          this is useful to be able to call C from Fortran
 *          and add the NEMO user interface
 *
 *	27-jun-91 recreated...		PJT
 *	24-may-92 added a 'c:' line as usage line (also conform MIRIAD)  PJT
 *      20-may-94 added hpux to not convert fortran external name
 *	 7-jan-00 F77_FUNC conversion		PJT
 *      30-mar-01 doc/typos	PJT
 *      20-jun-01 gcc3 
 */

#include <nemo.h>
#include <strlib.h>
#include <ctype.h>

string defv[] = {
    "in=???\n           Input source code to extract defv[] from",
    "out=-\n            Output nemo_main module [- means stdout]",
    "call=nemomain\n    Your fortran main subroutine (no args!)",
    "options=\n         define your own options (lower,upper,under)",
    "VERSION=2.0b\n     20-jun-01 PJT",
    NULL,
};

string usage="nemo_main() builder from (fortran) source code";
string e_usage=NULL;

#define MAXLINE 512

extern bool scanopt(string, string);

void mksymname(char *name, char *options)
{
    char *cp=name;
    int upper=0, lower=0, underend=0;

    if (name==0 || *name==0) return;

    if (options==0 || *options == 0) {
#if FORTRANIZE_LOWERCASE
	dprintf(0,"ftoc: using FORTRANIZE_LOWERCASE\n");
	lower=1;
#elif FORTRANIZE_UPPERCASE
	dprintf(0,"ftoc: using FORTRANIZE_UPPERCASE\n");
	upper=1;
#elif FORTRANIZE_LOWERCASE_UNDERSCORE
	dprintf(0,"ftoc: using FORTRANIZE_LOWERCASE_UNDERSCORE\n");
	lower=1; underend=1;
#elif FORTRANIZE_UPPERCASE_UNDERSCORE
	dprintf(0,"ftoc: using FORTRANIZE_UPPERCASE_UNDERSCORE\n");
	upper=1; underend=1;
#endif

    } else {
	dprintf(0,"%d %d %d\n",lower,upper,underend);
	if (scanopt(options,"lower")) lower=1;
	dprintf(0,"%d %d %d\n",lower,upper,underend);
	if (scanopt(options,"upper")) upper=1;
	dprintf(0,"%d %d %d\n",lower,upper,underend);
	if (scanopt(options,"under")) underend=1;
	dprintf(0,"%d %d %d\n",lower,upper,underend);
    }

    if (upper) {            /* convert to upper case */
        while (*cp) {
            if (islower(*cp))
                *cp = toupper(*cp);
            cp++;
        }
    }

    if (lower) {            /* convert to lower case */
        while (*cp) {
            if (isupper(*cp))
                *cp = tolower(*cp);
            cp++;
        }
    }

    if (underend) {         /* append underscore */
        while (*cp) 
            cp++;
        *cp++ = '_';
        *cp = '\0';
    }
}




void header(stream outstr)
{
    fprintf(outstr,"/* THIS FILE HAS BEEN CREATED BY %s - do not edit */\n\n",
            getargv0());
    fprintf(outstr,"#include <stdinc.h>\n\n");
    fprintf(outstr,"string defv[] = {\n");
}

void footer(stream outstr, string callname, string options)
{
    char name[MAXLINE];

    if (*callname)
        strcpy(name,callname);
    else
        strcpy(name,"nemomain");
    mksymname(name, options);
    
    fprintf(outstr,"    NULL,\n};\n\n");
    fprintf(outstr,"string usage=\"%s\";\n\n",
        (e_usage ? e_usage : "NEMO program with unknown intent"));
    fprintf(outstr,"extern void %s(void);\n",name);
    fprintf(outstr,"void nemo_main()\n{\n");
    fprintf(outstr,"    %s();\n}\n",name);
}

/* kept here for historic reasons, before usage of the FORTRAN.... macros */

void old_mksymname(char *name, char *machname)
{
    char *cp=name;
    int upper=0, underend=0;

    if (name==0 || *name==0) return;

    if (machname==0 || *machname == 0) {
#if defined(unicos) || defined(cray) || defined(cray2)
       upper=1;
#elif defined(hpux)
       ;
#else
       underend=1;
#endif    
    } else {
        if (streq(machname,"cray"))
            upper=1;
        else if (streq(machname,"vms") || streq(machname,"hpux"))
	    ;       /* do nothing */
	else if (streq(machname,"bsd"))
            underend=1;
        else {
            warning("MKSYSNAME: Machine type %s not uderstood, bsd assumed",
                    machname);
            underend=1;
        }

    }

    if (upper) {            /* convert to upper case */
        while (*cp) {
            if (islower(*cp))
                *cp = toupper(*cp);
            cp++;
        }
    }

    if (underend) {         /* append underscore */
        while (*cp) 
            cp++;
        *cp++ = '_';
        *cp = '\0';
    }
}

bool process(stream outstr,char *line)
{
    char *cp;
    static int doc_status = -1;     /* -1: no started 0:inside   1: done */
    static int lineno=0;

    lineno++;
    cp =line;

    if (*cp != 'c' && *cp != 'C') return TRUE;
    cp++;
    if (*cp == ':') {
        cp++;
        while (*cp == ' ' || *cp == '\t') cp++;     /* skipwhite */
        if (*cp==0 || *cp=='\n') return TRUE;
        if (cp[strlen(cp)-1] == '\n') /* patch possible newline from fgets() */
            cp[strlen(cp)-1] = '\0';
        e_usage = scopy(cp);
        return TRUE;
    }
    if (*cp == '+') {
        if (doc_status != -1) error("Illegal c+ at line %d\n",lineno);
        doc_status = 0;
        return TRUE;
    }
    if (*cp == '-') {
        if (doc_status != 0)  error("Illegal c- at line %d\n",lineno);
        return FALSE;
    }
    if (doc_status != 0) return TRUE;
    cp++;
    while (*cp == ' ' || *cp == '\t') cp++;  /* skipwhite */
    if (*cp == 0)                    /* skip 'blank' lines */
        return TRUE;
    if (cp[strlen(cp)-1] == '\n')       /* patch possible newline from fgets() */
        cp[strlen(cp)-1] = '\0';
    fprintf(outstr,"    \"%s\",\n",cp);
    return TRUE;
}

void nemo_main()
{
    stream instr, outstr;
    bool notdone=TRUE;
    char line[MAXLINE+1];

    instr = stropen(getparam("in"),"r");
    outstr = stropen(getparam("out"),"w");

    header(outstr);
    while(notdone) {
        if (fgets(line,MAXLINE,instr)==0)
            break;
        notdone = process(outstr,line);
    }
    if (notdone)
        error("Input file does not have properly closed defv section");
    strclose(instr);
    footer(outstr,getparam("call"),getparam("options"));
    strclose(outstr);
}
