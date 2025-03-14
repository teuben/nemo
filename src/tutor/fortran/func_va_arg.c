#include <stdio.h>
#include <stdarg.h>

int test_va_arg_f3_(int * s1, ...)
{
    va_list pa;
    char * myval;
    fprintf(stderr,">>>>>>>> s1 [%d]\n",*s1);
    va_start(pa, s1);
    myval = (char *)  va_arg(pa, char *);
    fprintf(stderr,">>>>>>>> myval [%s]\n",myval);
    return 1;
}
