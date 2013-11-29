/*
 * MKEPIDISK.C: set up a testdisk with particles in epicyclic motion
 *   28-aug-90  PJT  Created
 *              -kludge: some vrot= correction was kludged in from
 *               the start to make a vel.fie. look nicer
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>
#include <filestruct.h>

#include <snapshot/snapshot.h>
#include <snapshot/body.h>
#include <snapshot/put_snap.c>

string defv[] = {       /* DEFAULT INPUT PARAMETERS */
    "out=???\n            output file name (snapshot)",
    "rmin=0.0\n           inner disk radius ",
    "rmax=1.0\n           outer cutoff radius ",
    "nring=20\n           number of rings",
    "dr=0.05\n            Typical interstellar distance",
    "axirat=1.0\n         Axial ratio of resulting orbits",
    "omega=0.0\n          Relative pattern speed",
    "seed=0\n             usual random number seed ",
    "headline=\n          text headline for output ",
    "VERSION=1.1\n        PJT  5-nov-90",
    NULL,
};

real rmin, rmax, dr, omega, axirat;

int nbody, nring;

Body *btab;

main(argc, argv)
string argv[];
{
    initparam(argv, defv);

    rmin = getdparam("rmin");
    rmax = getdparam("rmax");
    dr = getdparam("dr");
    nring = getiparam("nring");
    omega = getdparam("omega");
    axirat = getdparam("axirat");

    set_xrandom(getiparam("seed"));
    testdisk();
    writegalaxy(getparam("out"), getparam("headline"));
}

/*
 * WRITEGALAXY: write galaxy model to output.
 */

writegalaxy(name, headline)
string name;
string headline;
{
    stream outstr;
    real tsnap = 0.0;
    int bits = MassBit | PhaseSpaceBit;

    if (! streq(headline, ""))
        set_headline(headline);
    outstr = stropen(name, "w");
    put_history(outstr);
    put_snap(outstr, &btab, &nbody, &tsnap, &bits);
    strclose(outstr);
}

testdisk()
{
    char *malloc();
    Body *bp;
    real dring, ccos, csin, cx, cy, r, phi, dphi, vrot;
    int i, k, n;
    double xrandom(), sin(), cos(), sqrt();

    dring = (rmax-rmin) / nring;
    for (r=rmin, i=0, nbody=0; i<nring; i++, r += dring) /* first pass to get nbody */
        nbody += (int)  (6.28*r/dr);
    dprintf(0,"Total of %d bodies in %d rings in the disk\n",nbody,nring);

    btab = (Body *) malloc(nbody * sizeof(Body));   /* allocate all particles */
    if (btab == NULL)
        error("mkepidisk: not enuf memory for %d particles\n",nbody);

    bp = btab;
    phi = 0;                    /* running angle */
    for (r=rmin, i=0; i<nring; i++, r += dring) {
        n = (int)  (6.28*r/dr);
        dphi = 6.28/n;          /* current step in angle for this ring */
        cx =  r * 2 / (1+axirat);
        cy =  r * 2 * axirat / (1+axirat);
        vrot = 4.0 /(1+4*r);     /* 'irrelevant' correction for diff.rot */
        for (k=0; k<n; k++, bp += 1) {
            Mass(bp) = 1.0 / nbody;           /* total mass = 1.0 */
            csin = sin(phi);
            ccos = cos(phi);
            phi += dphi;
            Pos(bp)[0] = cx * ccos;                /* set positions    */
            Pos(bp)[1] = cy * csin;
            Pos(bp)[2] = 0.0;
            Vel(bp)[0] = -cx * csin * vrot - omega*Pos(bp)[1];
            Vel(bp)[1] =  cy * ccos * vrot + omega*Pos(bp)[0];
            Vel(bp)[2] = 0.0;
        }
    }
}
