/*
 * f2c and c2f interface aids
 */

#include <stdinc.h> 
 
strcpy_f2c(s1,n1,s2,n2)
char *s1, *s2;
int n1, n2;
{
    int lnblnk=n2;

    if (s1==0 || s2==0) return;		/* valid strings */

    while(lnblnk > 0) {        /* count real string length */
        if (s2[lnblnk-1] != ' ') break;
        lnblnk--;
    }
    if (lnblnk <= 0) {         /* no (valid) string input */
        if (n1 > 0) *s1 = 0;
        return;
    }
    if (lnblnk >= n1) lnblnk = n1-1;
    strncpy(s1,s2,lnblnk);
    s1[lnblnk] = 0;
}

strcpy_c2f(s1,n1,s2,n2)
char *s1, *s2;
int n1, n2;
{
    int i, n = strlen(s2);

    if (n > n1) n = n1;
    strncpy(s1,s2,n);
    for (i=n; i<n1; i++) s1[i] = ' ';
}



