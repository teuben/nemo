/*
 * FORCE_MATH: this curious function is never called.  However, its mere
 * existence is sufficient to fool ld into loading the math functions.
 * Useful to include by code which uses the dynamic object loader.
 *
 */

static force_math()
{
    real cbrt(), sqrt();
    real sin(), cos(), asin(), acos();
    real tan(), atan(), atan2();
    real exp(), dex(), log(), log10(), pow();
    real fabs(), floor(), ceil(), rint();

    int spline();
    real log2(), erf(), sqr();

    (void) cbrt(1.0);
    (void) sqrt(1.0);
    (void) sin(1.0);
    (void) cos(1.0);
    (void) asin(1.0);
    (void) acos(1.0);
    (void) tan(1.0);
    (void) atan(1.0);
    (void) atan2(1.0, 1.0);
    (void) exp(1.0);
    (void) dex(1.0);
    (void) log(1.0);
    (void) log10(1.0);
    (void) log2(1.0);
    (void) pow(1.0, 1.0);
    (void) fabs(1.0);
    (void) floor(1.0);
    (void) ceil(1.0);
    (void) rint(1.0);

    (void) spline();	/* not the right arguments though */
    (void) sqr(1.0);
    (void) erf(1.0);
}
