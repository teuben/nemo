#include <stdio.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "vtc.h"
#include "vtclocal.h"

static double t0;

#if 0

/* CPU time */
void
vtc_init_cputime(void)
{
    struct rusage x;
    double sec,microsec;

    getrusage(RUSAGE_SELF,&x);
    sec = x.ru_utime.tv_sec + x.ru_stime.tv_sec ;
    microsec = x.ru_utime.tv_usec + x.ru_stime.tv_usec ;
    t0 = (sec + microsec / 1000000.0);
}

/* CPU time - t0 */
double
vtc_get_cputime(void)
{
    struct rusage x;
    double sec,microsec;

    getrusage(RUSAGE_SELF,&x);
    sec = x.ru_utime.tv_sec + x.ru_stime.tv_sec ;
    microsec = x.ru_utime.tv_usec + x.ru_stime.tv_usec ;
    return (sec + microsec / 1000000.0 - t0);
}

#else

/* wallclock time */
void
vtc_init_cputime(void)
{
    struct timeval x;

    gettimeofday(&x, NULL);
    t0 = (x.tv_sec + x.tv_usec / 1000000.0);
}

/* wallclock time - t0 */
double
vtc_get_cputime(void)
{
    struct timeval x;

    gettimeofday(&x, NULL);
    return(x.tv_sec + x.tv_usec / 1000000.0 - t0);
}

#endif

static int cputime_on = 0;

void
vtc_turnon_cputime(void)
{
    cputime_on = 1;
}

void
vtc_turnoff_cputime(void)
{
    cputime_on = 0;
}

void
vtc_print_cputime(char *comment)
{
    if (cputime_on) {
	fprintf(stderr, "%s %6.5f\n", comment, vtc_get_cputime());
    }
}
