/*
 *  select: table row selection routines
 *
 *	18-dec-99	NULL->0			pjt
 *	20-jun-01 	gcc3			pjt
 *
 */
#include <stdinc.h>
#include <getparam.h>

   /* allowed numeric comparisons 'a OPER b' with the following OPER's */

#define ISBAD   0      
#define ISLE    1   /* <= */
#define ISLT    2   /* <  */
#define ISGE    3   /* >= */
#define ISGT    4   /* >  */
#define ISEQ    5   /* == */
#define ISNE    6   /* != */

local int substitute(string s, string *cnam, string *cval, 
	real *a, real *b);


extern void strinsert(char *a, char *b, int n);

bool num_select(string *cnam, string *cval, string *sels)
{
    string *s;
    int retval = 1;     /* set return value initially to TRUE */
    real a,b;

    if (sels==NULL || *sels==NULL) return TRUE;

    for (s=sels; *s; s++) {
        switch(substitute(*s, cnam, cval, &a, &b)) {
            case ISLT:  retval &= (a < b);      
                        break;
            case ISLE:  retval &= (a <= b);      
                        break;
            case ISGT:  retval &= (a > b);      
                        break;
            case ISGE:  retval &= (a >= b);      
                        break;
            case ISEQ:  retval &= (a == b);      
                        break;
            case ISNE:  retval &= (a != b);      
                        break;
            default:    error("Illegal opcode in %s",*s);
        }
        if (retval==0) return FALSE;        /* bail out if FALSE already */
    }
    return TRUE;        /* if gotten here, all is TRUE and return as such */
}

local int substitute(string s, string *cnam, string *cval, real *a, real *b)
{
    char sel[128], var[128];
    char *vp, *cp = sel;
    string *cn, *cv;
    int vlen, retval;
    bool found;
    
    strcpy(sel,s);  /* patching is done in a local variable */
    for (;;) {      /* loop over selections to substitute variables */
        cp = strchr(cp,'%');        /* search for reference char */
        if (cp==NULL) break;        /* done when no more found */
        strcpy(var,cp+1);           /* copy remainder in var */
        vp = var;
        while (*vp != 0 && *vp != ' ')   /* FIX THIS !!!!! */
            vp++;
        *vp = 0;
        vlen = strlen(var) + 1;     /* count length of '%VARNAME' */

        found = FALSE;        
        for (cn=cnam, cv=cval; *cn; cn++, cv++) {
            if (streq(var,*cn)) {
                found = TRUE;
                strinsert(cp,*cv,vlen);
                break;
            }
        }
        if (!found) error("No known variable %%%s",var);
        
    }

    if ((cp=strchr(sel,'<'))) {               /* find what the operand is */
        *cp++ = 0;
        if (*cp == '=') {
            retval = ISLE;
            *cp++ = 0;
        } else
            retval = ISLT;
    } else if ((cp=strchr(sel,'>'))) {
        *cp++ = 0;
        if (*cp == '=') {
            retval = ISGE;
            *cp++ = 0;
        } else
            retval = ISGT;
    } else if ((cp=strchr(sel,'='))) {
        *cp++ = 0;
        if (*cp == '=') {
            retval = ISEQ;
            *cp++ = 0;
        } else
            retval = ISBAD;
    } else if ((cp=strchr(sel,'!'))) {
        *cp++ = 0;
        if (*cp == '=') {
            retval = ISNE;
            *cp++ = 0;
        } else
            retval = ISBAD;
    } else
        retval = ISBAD;

    if (retval==ISBAD) return ISBAD;

    if (nemoinpr(sel,a,1) != 1) warning("Error parsing %s",sel);
    if (nemoinpr(cp, b,1) != 1) warning("Error parsing %s",cp);
        
    return retval;
}
