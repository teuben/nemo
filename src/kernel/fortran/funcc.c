#include <nemo.h>

#define funcc F77_FUNC(funcc,FUNCC)
#define funcf F77_FUNC(funcf,FUNCF)
        
extern double funcf(double *);
   
double funcc(double * xp) {
    double x = *xp;
    dprintf(1,"FuncC: x=%g\n",x);
    if (x < 1.0)
        return  sqrt(x);
    else {
        x = 1.0/x;
        return  funcf(&x);
    }
}



