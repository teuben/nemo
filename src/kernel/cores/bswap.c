/* 
 *   BSWAP: swap bytes - called by lower level filestruct routines
 *          if a byte swapped data set is detected at input.
 *		dat	pointer to the data
 *		len	item length in bytes, this amount will be swapped
 *		cnt	total number of data-items of length 'len' to be
 *			swapped
 * 
 *   TESTBED: can also be used to tweek this routine and measure 
 *            simple CPU operations.
 *   BUGS: cannot check if dat is long enough (len*cnt)
 *   NOTES:  An IBM RISC  is about 40% slower than a SPARC1
 *          
 *   History:
 *	12-oct-90  V1.0 PJT	Written
 *	20-nov-91  V1.1 PJT  	NEMO V2.x nemo_main interface, added out=
 *	21-nov-91          	TOOLBOX version called bswap for Ultrix-SUN 
 *                      word swap for Stephen White; 'dd conv=swab' will do
 *      27-nov-91  fixed typo (missing declared k; l->k)
 *	25-feb-92  happy gcc2.0
 *	20-nov-92  TOOLBOX, not TESTBED
 *	25-feb-94  argument now (void *) from (char *) for ansi
 *       9-apr-94  also handle pipes in TOOLBOX section
 *       4-aug-94  warn about small files where t1==t0 etc. -- TOOLBOX only
 *	23-feb-97  allow offset until which no swapping done
 *	15-dec-98  TOOLBOX : allow streaming mode
 *      12-sep-01  file_size
 */

#include <stdinc.h>

void bswap(void *vdat, int len, int cnt)
{
    char tmp, *dat = (char *) vdat;
    int k;

    if (len==1)
	return;
    else if (len==2)
        while (cnt--) {
            tmp = dat[0];  dat[0] = dat[1];  dat[1] = tmp;
            dat += 2;
        }
    else if (len==4)
        while (cnt--) {
            tmp = dat[0];  dat[0] = dat[3];  dat[3] = tmp;
            tmp = dat[1];  dat[1] = dat[2];  dat[2] = tmp;
            dat += 4;
        }
    else if (len==8)
        while (cnt--) {
            tmp = dat[0];  dat[0] = dat[7];  dat[7] = tmp;
            tmp = dat[1];  dat[1] = dat[6];  dat[6] = tmp;
            tmp = dat[2];  dat[2] = dat[5];  dat[5] = tmp;
            tmp = dat[3];  dat[3] = dat[4];  dat[4] = tmp;
            dat += 8;
        }
    else {  /* the general SLOOOOOOOOOWE case */
        for(k=0; k<len/2; k++) {
            tmp = dat[k];
            dat[k] = dat[len-1-k];
            dat[len-1-k] = tmp;
        }
    }
}

#if defined(TOOLBOX)

#include <getparam.h>

string defv[] = {
    "in=???\n       File to read and swap bytes of",
    "out=\n         Filename if the swapped file to be output",
    "len=2\n        Itemlength in bytes during swapping",
    "oneswap=t\n    One swap call? (or many) - for testing only", 
    "offset=0\n     Offset (in bytes) before which no swapping done",
    "VERSION=1.4\n  15-dec-98 PJT",
    NULL,
};

string usage="swap bytes in a file";

extern int file_size(string);

void nemo_main(void)
{
    stream instr, outstr;
    string fname;
    char *data;
    real t0, t1, t2, rspeed=0, wspeed=0;
    int i, len, cnt, offset;
    bool onetrip;

    fname = getparam("in");
    len = getiparam("len");
    offset = getiparam("offset");
    onetrip = getbparam("oneswap");
    if (hasvalue("out"))
        outstr = stropen(getparam("out"),"w");
    else {
        warning("No swapped output file created");
        outstr = NULL;
    }

    instr = stropen(fname,"r");
    cnt = file_size(fname);         /* size of the file in bytes */
 if (cnt < 0) {          /* streaming mode */
    dprintf(0,"Streaming mode\n");
    data = (char *) allocate(len);
    while (fread(data,len,1,instr) == 1) {
        bswap(data,len,1);
        fwrite(data,1,len,outstr);
    }
   strclose(outstr);
 } else {
    if (cnt % len)
        warning("Filesize is not a multiple of itemlength: %d / %d", cnt, len);
    cnt /= len;                     /* cnt is now number of items of 'len' */
    if (cnt==0) error("File %s too small; file_size() returned 0",fname);
    data = (char *) allocate(cnt*len);
    t0 = 60 * cputime();
    if (cnt != fread(data,len,cnt,instr)) error("Error reading %s",fname);
    t1 = 60 * cputime();
    if (onetrip)
        bswap(data,len,cnt);
    else {
        i = cnt;
        while (i--)
        bswap(data,len,1);
    }
    t2 = 60 * cputime();
    if (t1==t0) warning("Inaccurate measurement for reading");
    else rspeed=cnt*len/(t1-t0)/1048576.0;
    if (t1==t2) warning("Inaccurate measurement for writing");
    else wspeed=cnt*len/(t2-t1)/1048576.0;
    dprintf(1,"Read %d: %f s; swap %d * %d bytes: %f s; %g %g Mswap\n",
       		cnt*len, t1-t0, cnt, len, t2-t1, 
		rspeed,wspeed);
    printf("%d %d %d %g %g\n",
                cnt,len,(onetrip?1:0),
		rspeed,wspeed);
    if (outstr) {
        if (cnt != fwrite(data,len,cnt,outstr)) 
            error("Error writing %s",getparam("out"));
        strclose(outstr);
    }
  }
}
#endif

