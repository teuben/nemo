#include <nemo.h>

string defv[] = {
    "x=4.0\n        number to test",
    "VERSION=1.0\n  7-jan-00  PJT",
    NULL,
};

string usage = "F77_FUNC macro tester";



#define funcc F77_FUNC(funcc,FUNCC)
#define funcf F77_FUNC(funcf,FUNCF)
        
extern double funcf(double *);
extern double funcc(double *);


nemo_main()
{
    double x = getdparam("x");

    printf("funcc(x) = %g\n",funcc(&x));
    printf("funcf(x) = %g\n",funcf(&x));
}


