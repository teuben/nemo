/*
 * log.c:  logarithmic triaxial potential
 *	dec-93    allowed r_c=0 exception, by setting v_0^2 = 2*m_c
 */

/*CTEX
 *    {\bf potname=log
 *       potpars={\it $\Omega,M_c,r_c,q$}} 
 *
 * Logarithmic Potential (BT, pp.45, eq. 2.54 and eq. 3.77
 *
 * $$
 *    \Phi = {1\over 2} v_0^2
 *                     \ln{ \left( r_c^2 + r^2 \right) }
 * $$
 *
 * with $ M_c \equiv {1\over 2} r_c v_0^2 $ defined as the ``core mass''.
 */
                      
 
/* default parameters */

static double omega = 0.0;	/* pattern speed */
static double mc = 1.0;		/* core mass */
static double rc = 1.0;		/* core radius */
static double q  = 1.0;		/* asphericity (q<=1)	*/
static double r;		/* default: r=q (prolate) */

static double mor, r2, iq2[3];  /* scratch variables  */

extern double sqr(), log();

void inipotential (npar, par, name)
int    *npar;
double par[];
char *name;
{
    int i;

    if (*npar>0) omega = par[0];
    if (*npar>1) mc = par[1];
    if (*npar>2) rc = par[2];
    if (*npar>3) q = par[3];
    if (*npar>4)
	r = par[4];		/* if it's there: make triaxial */
    else
    	r = q;			/* default: prolate */

    dprintf (1,"INI_POTENTIAL Logarithmic potential (prolate) %s\n",name);
    dprintf (1,"  Parameters : Pattern Speed = %f\n",omega);
    dprintf (1,"  m_c, r_c, q= %f %f %f %f\n",mc,rc,q,r);

    if (rc>0)    
        mor = mc/rc;
    else
    	mor = mc;
    r2 = rc*rc;
    iq2[0] = 1.0;
    iq2[1] = 1.0/sqr(q);
    iq2[2] = 1.0/sqr(r);
    par[0] = omega;
}
    

void potential (ndim,pos,acc,pot,time)
int    *ndim;
double pos[], acc[], *pot, *time;
{
    double rad, f;
    int    i;
	
    for (i=0, rad=r2; i<*ndim; i++)
        rad += sqr(pos[i])*iq2[i];
    *pot = mor * log(rad);
    f = -2.0*mor/rad;
    for (i=0; i<*ndim; i++)
        acc[i] = f*pos[i]*iq2[i];
}
