/*
 * UNFIO:  C routines that allow access to binary files created
 *	   with fortran in UNFORMATTED mode
 *         Only works with UNIX f77-like compilers, don't even think
 *         of trying this on Cray, VMS or Microsoft MSDOS or any non-IEEE
 *         machines (:-) Well, it's just a matter of #ifdeffing this....
 *
 *	may-94	created in a hurry for rvcsnap(1NEMO)
 *   20-may-94  fixed bug in computing  'n' in display(); added warning EOF.
 *   22-may-95  added some comments, extra check at end of unfread
 *    8-aug-95  allow swapped reads, and added some dprintf()'s
 *   21-jun-97  added select= keyword (1-based)
 *    7-feb-98  2.0  added out= to replace the 'cf' utility (src/image/fits/)
 *   19-mar-99  2.1  out= now also listens to select= (albeit slowly)
 *   20-jun-01  gcc3
 *   
 */

#include <stdinc.h>
#include <unfio.h>

static bool do_swap;    /* (re)set in a call to unfswap() */

extern void bswap(void *vdat, int len, int cnt);


/*
 * unfswap: set swapping mode, by default, no swapping performed
 *
 */

int unfswap(bool swap)
{
    do_swap = swap;
    return 0;
}

/*
 * unfscan: read another block; return the number of bytes in that block
 *          or 0 if EOF or another error (e.g. if the file doesn't seem
 *          to look like a blocked file
 */

int unfscan(stream fp)
{
    int n, size0, size1;
    
    n = fread(&size0,sizeof(int),1,fp);
    if (n!=1) return 0;
    if (do_swap) bswap(&size0, sizeof(int), 1);
    dprintf(2,"unfread: header %d\n",size0);

    fseek(fp,size0,1);

    n = fread(&size1,sizeof(int),1,fp);
    if (n!=1) return -1;
    if (do_swap) bswap(&size1, sizeof(int), 1);
    dprintf(2,"unfread: trailer %d\n",size1);
    if (size0 != size1) return -2;

    return size0;

}

/*
 * unfread: read another block; return the number of bytes in that block
 *          or 0 if EOF or another error (e.g. if the file doesn't seem
 *          to look like a blocked file. Return the data in the buffer
 *          pointed to by 'buf'. Fatal error if buffer not big enough.
 */

int unfread(stream fp, char *buf, int bufsize)
{
    int n, size, size1;

    n = fread(&size,sizeof(int),1,fp);
    if (n!=1) return 0;
    if (do_swap) bswap(&size, sizeof(int), 1);
    dprintf(2,"unfread: header %d\n",size);

    if (size > bufsize) {
        error("unfread: buffersize %d too small to read all %d data",
                size,bufsize);
	return 0;
    } else
        dprintf(2,"unfread: header %d\n",size);

    n = fread(buf, sizeof(char), size, fp);
    dprintf(2,"unfread: data %d\n",n);
    if (n != size) return 0;

    n = fread(&size1, sizeof(int), 1, fp);
    if (n != 1) return 0;
    if (do_swap) bswap(&size1, sizeof(int), 1);
    dprintf(2,"unfread: trailer %d\n",size1);
    if (size1 != size) 
        warning("Reading block size_s=%d start_e=%d",size,size1);

    return size;
}


#ifdef TOOLBOX

#include <getparam.h>

string defv[] = {
        "in=???\n           input",
        "out=\n             output",
        "block=\n           which block to display (0=all, default: scan)",
        "type=i\n           format {int, float, double}",
        "select=\n          which elements from this block to select (1=first)",
        "format=\n          printf format to use if block display [%d,%g]",
        "count=f\n          display element counter too?",
        "maxbuf=10000\n     buffersize in bytes, to read a block",
	"swap=f\n           swapped read?",
        "VERSION=2.1a\n	    21-jun-01 PJT",
        NULL,
};

string usage = "access fortran unformatted I/O files";

void my_display(int, char *, string, string, int, int *, bool, bool, stream);

void nemo_main()
{
    int n, iblock, block=0;
    stream instr = stropen(getparam("in"),"r");
    stream ostr = NULL;
    bool Qdisp = hasvalue("block");
    bool Qcount = getbparam("count");
    bool Qswap = getbparam("swap");
    bool Qsel = hasvalue("select");
    string type = getparam("type");
    string fmt = getparam("format");
    int maxbuf = getiparam("maxbuf");
    char *buf = (char *) allocate(maxbuf);
    int nidx, *idx = NULL;

    if (Qdisp) iblock = getiparam("block");
    if (hasvalue("out")) ostr = stropen(getparam("out"),"w");
    unfswap(Qswap);
    if (Qsel) {
        nidx = maxbuf/sizeof(int);
        idx = (int *) allocate(maxbuf);
        nidx = nemoinpi(getparam("select"),idx,nidx);
    } else
        nidx = 0;

    if (Qdisp) {
        for (;;) {
            block++;
            if (iblock==0 || iblock==block) {
                n = unfread(instr,buf,maxbuf);
                if (n<1) break;
                my_display(n, buf, type, fmt, nidx, idx, Qcount, Qswap, ostr);
            } else {
                n = unfscan(instr);
                if (n<1) {
		  if (iblock>=block) 
		     warning("Requested block (%d) too large (%d)",iblock,block);
		  break;
		}
            }
        }
    } else {
        for (;;) {
            n = unfscan(instr);
            if (n<1) break;
            block++;
            dprintf(0,"%d   %d\n",block,n);
        }
    }
}

void my_display(int n, char *buf, string type, string fmt, 
                int nidx, int *idx, bool count, bool swap, stream ostr)
{
    double *dp;
    float *fp;
    int *ip;
    int i, j, nmax;
    char format[32];
    bool Qdef, Qidx;

    Qdef = (fmt==NULL || *fmt==0);
    Qidx = (nidx > 0);

    if (ostr && !Qidx) {
        if (fwrite(buf,1,n,ostr) != n)
            error("Error writing block");
        return;
    }


    switch (*type) {
    case 'i':
        ip = (int *) buf;
	n /= sizeof(int);
	if (swap) bswap(buf, sizeof(int), n);
	strcpy(format, Qdef ? "%d" : fmt);
        nmax = Qidx ? nidx : n;
        for (i=0; i<nmax; i++) {
            if (Qidx) {
                j = idx[i]-1;
                if (j>=n || j<0) continue;
            } else
                j = i;
            if (ostr)
                fwrite(&ip[j],1,sizeof(int),ostr);
            else {                
                if (count) printf("%d ",j+1);
	        printf(format,ip[j]);
	        printf("\n");
	    }
        }
        break;
    case 'f':
        fp = (float *) buf;
	n /= sizeof(float);
	if (swap) bswap(buf, sizeof(float), n);
	strcpy(format, Qdef ? "%g" : fmt);
        nmax = Qidx ? nidx : n;
        for (i=0; i<nmax; i++) {
            if (Qidx) {
                j = idx[i]-1;
                if (j>=n || j<0) continue;
            } else
                j = i;
            if (ostr)
                fwrite(&fp[j],1,sizeof(int),ostr);
            else {
                if (count) printf("%d ",j+1);
                printf(format,fp[j]);
                printf("\n");
            }
        }
        break;
    case 'd':
        dp = (double *) buf;
	n /= sizeof(double);
	if (swap) bswap(buf, sizeof(double), n);
	strcpy(format, Qdef ? "%g" : fmt);
        nmax = Qidx ? nidx : n;
        for (i=0; i<nmax; i++) {
            if (Qidx) {
                j = idx[i]-1;
                if (j>=n || j<0) continue;
            } else
                j = i;
            if (ostr)
                fwrite(&dp[j],1,sizeof(int),ostr);
            else {
                if (count) printf("%d ",j+1);
                printf(format,dp[j]);
                printf("\n");
            }
        }
        break;
    default:
        warning("Type %s not supported (i,f,d)",type);
    }
}
#endif
