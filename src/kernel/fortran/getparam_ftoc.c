/* ? f to c only has access to the getparam, since initparam etc. is
 *   done in C stuff 
 *
 *   BSD interface test	
 *	26-jun-91  quick hack together - i know, it's bad...	PJT
 *      22-feb-94  ANSI headers
 *	14-apr-95  no more ARGS
 *	 6-jan-00  converted to F77_FUNC macros
 */

#include <nemo.h>
#include <getparam.h>

local int len1(char *, int);
local int stof(string, string, int);
local char *stoc(string, int);

/* FORTRAN interfaces defined here:
 *    cvalue = getparam(key)
 *    ivalue = getiparam(key)
 *    dvalue = getdparam(key)
 */

#define fgetparam  F77_FUNC(getparam, GETPARAM)
#define fgetdparam F77_FUNC(getdparam,GETDPARAM)
#define fgetiparam F77_FUNC(getiparam,GETIPARAM)

void fgetparam(retval,retval_len, key,key_len)
    string key,     retval;
    int    key_len, retval_len;
{
    char *key_c, *val_c;

    key_c = stoc(key, key_len); 
    val_c = getparam(key_c);
    if(stof(val_c,retval,retval_len))
        warning("GETPARAM: Could not store full keyword value for %s:\n   %s\n",
                key_c, val_c);
    free(key_c);
    free(val_c);
}    

double fgetdparam(key, key_len)
    string key;
    int    key_len;
{
    double val;
    char *key_c;

    key_c = stoc(key, key_len);
    val = getdparam(key_c);
    free(key_c);
    return val;
}

int fgetiparam(key, key_len)
    string key;
    int    key_len;
{
    int val;
    char *key_c;

    key_c = stoc(key, key_len);
    val = getiparam(key_c);
    free(key_c);
    return val;
}


local string stoc(s, s_len)	/* take fortran string, return C string */
    string s;
    int s_len;
{
    char *cp;
    int l;


    l = len1(s,s_len);
    cp = allocate(l+1);     /* allocate 1 extra for the NULL */
    strncpy(cp,s,l);
    cp[l] = '\0';           /* and patch the NULL */
    dprintf(1,"STOC: s_len=%d -> l=%d s=%s\n",s_len,l,cp);
    return(cp);
}

				/* convert C string to fortran string */
local int stof(s, fs, fs_len)	/* return 0 if OK, 1 if not s not fitted in fs */
    char *s;            /* input C string */
    char *fs;           /* pointer to FORTRAN string (overwritten) */
    int fs_len;         /* declared length of CHARACTER variable */
{
    int i, l, retval=0;

    l = strlen(s);
    dprintf(1,"STOF: l=%d -> fs_len=%d s=%s\n,",l,fs_len,s);
    if (l>fs_len) {
        retval=l;
        l = fs_len;
    }
    for (i=0; i<l; i++)         /* copy part before the NULL */
        *fs++ = *s++;
    while(i++<fs_len)		/* blank fill the full fortran part */
        *fs++ = ' ';    
    if (retval) warning("STOF: %d/%d: could only save part of your string",
                            retval,fs_len);
    return(retval);
}

local int len1(char *s, int s_len)   /* true used length of a fortran character string */
{
    char *cp;
    int l = s_len;

    cp = &s[s_len-1];       /* point to last element */
    while (*cp == ' ' && l>0) {
        l--;
        cp--;
    }
    return l;
}
    

