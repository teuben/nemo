/*
 * C counterpart of nbody_io.f
 *
 * Assuming this routine emulated unformatted fortran I/O properly
 * (works on most BSD based fortran implementations, but probably
 *  not on odd machines like VMS, Cray etc.)
 * this is the only way how to do I/O on swapped binary files
 * and INTEGER*2 datasets in 'older' data.
 *
 *  8-aug-95    created             pjt
 * 20-jun-01    gcc3 casting        pjt
 */

#include <stdinc.h>
#include <unfio.h>


static stream unit3 = NULL;
static int sizeof_name = sizeof(int);
static int sizeof_data = sizeof(float);
static int Qswap;

extern void bswap(void *vdat, int len, int cnt);

void nb3open_c(string fname, int ilen, bool swap)
{

    if (unit3) error("unit3 already open - cannot handle two");

    unit3 = stropen(fname,"r");
    sizeof_name = ilen;
    Qswap = swap;
    unfswap(Qswap);
    dprintf(0,"nb3header_c: Opening %s using INTEGER*%d in %s mode\n",
            fname, sizeof_name, Qswap ? "swapped" : "native");
}

void nb3header_c(int *n, int *model, int *nrun, int *nk)
{
    int nread, header[4];

    unfswap(Qswap);
    if (*nk == 0) {
        /* READ (3, ERR=99, END=99)  n, model, nrun */
        nread = unfread(unit3, (char *)header, 3*sizeof(int));
	if (nread < 1) {
	    *n = 0;
	    return;
	}
        if (Qswap) bswap(header, sizeof(int), 3);
    } else {
        /* READ (3, ERR=99, END=99)  n, model, nrun, nk */
        nread = unfread(unit3, (char *)header, 4*sizeof(int));
	if (nread < 1) {
	    *n = 0;
	    return;
	}
        if (Qswap) bswap(header, sizeof(int), 4);
        *nk = header[3];
    }
    *n = header[0];
    *model = header[1];
    *nrun = header[2];
    return;
}

void nb3data_c(int *n, int *nk, float *a, 
              float *body, float *xs, float *xdot, int *name)
{
    int i, j, k, nread, nbuf, n1, n2, count;
    char *buf;
    float *fbuf;
    short *sbuf;
    int   *ibuf;

    /*
     * READ (3)  (a(k),k=1,nk), (body(j),j=1,n),
     *           ((xs(k,j),k=1,NDIM),j=1,n), ((xdot(k,j),k=1,NDIM),j=1,n),
     *           (name(j),j=1,n)
     */

    if (*n <= 0 || *nk <= 0) return;
    unfswap(Qswap);
    n1 = ((*nk) + 7 * (*n));  /* size of all REAL arrays */
    n2 = (*n) ;               /* size of the trailing INTEGER*? array */
    nbuf = n1 * sizeof_data + n2 * sizeof_name ;
    buf = (char *) allocate(nbuf);
    nread = unfread(unit3,buf,nbuf);
    if (Qswap) {		/* swap data if needed */ 
        bswap(&buf[0],                sizeof_data, n1);
        bswap(&buf[n1 * sizeof_data], sizeof_name, n2);
    }
    /* transfer all float's */
    fbuf = (float *) &buf[0];
    for(count=0, k=0; k < (*nk); k++)       /* A header array */
        a[k] = fbuf[count++];
    for(j=0; j < (*n); j++)                 /* masses */
        body[j] = fbuf[count++];
    for(i=0, j=0; j < (*n); j++)            /* positions */
        for (k=0; k < 3; k++)
            xs[i++] = fbuf[count++];
    for(i=0, j=0; j < (*n); j++)            /* velocities */
        for (k=0; k < 3; k++)
            xdot[i++] = fbuf[count++];
    if (sizeof_name == 2) {                 /* INTEGER*2 name */
        sbuf = (short *) &fbuf[count];
        for(count=0, j=0; j < (*n); j++)
            name[j] = sbuf[count++];
    } else {                                /* INTEGER*4 name */
        ibuf = (int *) &fbuf[count];
        for(count=0, j=0; j < (*n); j++)
            name[j] = ibuf[count++];
    }
    free(buf);
}

void nb3close_c()
{
    if (unit3==NULL) {
        warning("unit3 was already closed");
    } else {
        strclose(unit3);
        unit3 = NULL;
    }
}
