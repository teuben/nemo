/*
 *	24-oct-90	Created for new 'NEMO'		PJT
 */

#include <stdinc.h>
#include <getparam.h>
#include <filestruct.h>

string defv[] = {
    "n=\n    	   Test number",
    "VERSION=1.0\n 24-oct-90 PJT",
    NULL,
};

nemo_main()
{
    string name;
    char fname[32];
    int n,i;
    stream str;

    for(;;) {
        promptparam("n","Enter value for n: ");
        n = getiparam("n");
        if (n<0) break;
        printf("N=%d\n",n);
    }
}
