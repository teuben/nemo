/*  f2c_string()
 *  
 * POTENTIAL F_TO_C
 *
 * stub in order for Fortan to call our C routines.
 *
 */


#include <stdio.h>

extern char *malloc();

char *f2c_string();

void inipotential_(int *npar, double *par, char * name, int namelen)
{
    char *local_name = f2c_string(name, namelen);
    
    inipotential(npar, par, local_name);
}

void potential_ (int ndim,double *pos,double *acc,double *pot,double *time)
{
    potential(ndim,pos,acc,pot,time);
}


char *f2c_string(char *s, int slen)
{
    char *cp;
    int i;

    if (s==NULL || *s==0) return NULL;
    if (slen<0) return NULL;

    for (i=slen-1; i>=0; i--)
        if (s[i] != ' ') break;
    i++;

    cp = (char *) malloc(i+1);
    strncpy(cp, s, i);
    cp[i] = 0;
    return cp;
}
