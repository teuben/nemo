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
 *   See also:  RFC 1832 e.g.: http://www.faqs.org/rfcs/rfc1832.html
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
 *      30-sep-03  testing memcpy, and improved the testing
 *      20-sep-05  little and big endian versions
 *      14-may-12  optionally use the ffswapX routines from cfitsio
 */

//#define HAVE_CFITSIO
//#define HAVE_FFSWAP

#include <stdinc.h>
#if defined(HAVE_CFITSIO)
#include "fitsio2.h"
#endif

void bswap(void *vdat, int len, int cnt)
{
    char tmp, *dat = (char *) vdat;
    int k;
#if defined(HAVE_FFSWAP)
    if (len==1)
	return;
    else if (len==2)
      ffswap2((short int *)dat,cnt);
    else if (len==4)
      ffswap4((int *) dat,cnt);
    else if (len==8)
      ffswap8((double *)dat,cnt);
    else {  /* the general SLOOOOOOOOOWE case; should never happen */
        for(k=0; k<len/2; k++) {
            tmp = dat[k];
            dat[k] = dat[len-1-k];
            dat[len-1-k] = tmp;
        }
    }
#else
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
#endif
}

/*
 * bswap_bigend:   bswap only if the source data was big endian
 */

void bswap_bigend(void *vdat, int len, int cnt)
{
#if !defined(WORDS_BIGENDIAN)
  bswap(vdat,len,cnt);
#endif
}

/*
 * bswap_litend:   bswap only if the source data was little endian
 */

void bswap_litend(void *vdat, int len, int cnt)
{
#if defined(WORDS_BIGENDIAN)
  bswap(vdat,len,cnt);
#endif
}

#if defined(TOOLBOX)

#include <getparam.h>

string defv[] = {
    "in=???\n       File to read and swap bytes of",
    "out=\n         Filename if the swapped file to be output",
    "len=2\n        Itemlength in bytes during swapping",
    "oneswap=t\n    One swap call? (or many) - for testing only", 
    "offset=0\n     Offset (in bytes) before which no swapping done",
    "endian=\n      assume 'Little' (l) or 'Big' (b) endian input file",
    "memcpy=f\n     Testing swapping double another way with memcpy",
    "repeat=0\n     How many times to repeat the swapping (speed testing)",
    "VERSION=2.0\n  14-may-2012 PJT",
    NULL,
};

#if defined(HAVE_FFSWAP)
string usage="swap bytes in a file - fast";
#else
string usage="swap bytes in a file - generic";
#endif

string cvsid="$Id$";

extern int nemo_file_size(string);


void byteswap_doubles(double *a)
{
  unsigned char b[8],c[8];
  memcpy(b,a,8); 
  c[0]=b[7]; /* swap data around */
  c[1]=b[6];
  c[2]=b[5];
  c[3]=b[4];
  c[4]=b[3];
  c[5]=b[2];
  c[6]=b[1];
  c[7]=b[0];
  memcpy(a,c,8);
}

typedef void  (*bswap_proc)(void *, int, int);


void nemo_main(void)
{
    stream instr, outstr;
    string fname, endian;
    char *data, *dp;
    real t0, t1, t2, rspeed=0, wspeed=0;
    int i, len, cnt, offset, repeat, nrepeat;
    bool onetrip;
    bool Qmemcpy = getbparam("memcpy");
    bswap_proc bptr;

#if defined(HAVE_FFSWAP)
    dprintf(1,"FFSWAP enabled\n");
#endif

    fname = getparam("in");
    len = getiparam("len");
    offset = getiparam("offset");
    nrepeat = repeat = getiparam("repeat");
    onetrip = getbparam("oneswap");
    if (hasvalue("endian")) {
      endian = getparam("endian");
      if (*endian == 'l' || *endian == 'L')
	bptr = bswap_litend;
      else if (*endian == 'b' || *endian == 'B')
	bptr = bswap_bigend;
      else
	error("Bad endian=%s; need 'big' or 'little' endian",endian);
    } else
      bptr = bswap;
    if (hasvalue("out"))
        outstr = stropen(getparam("out"),"w");
    else {
        warning("No swapped output file created");
        outstr = NULL;
    }
    if (Qmemcpy)
      if (len != 8) {
	warning("Cannot do memcpy, len=%d, needs to be 8",len);
	Qmemcpy = FALSE;
      } else
	warning("memcpy mode");

    instr = stropen(fname,"r");
    cnt = nemo_file_size(fname);         /* size of the file in bytes */

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
    if (cnt==0) error("File %s too small; nemo_file_size() returned 0",fname);
    data = (char *) allocate(cnt*len);
    t0 = 60 * cputime();
    if (cnt != fread(data,len,cnt,instr)) error("Error reading %s",fname);
    t1 = 60 * cputime();
    do {
      if (onetrip) {
	/* by doing a one-trip, speed is about 160 ; 275 for double */
        bswap(data,len,cnt);
      } else {
	/*
	 * by incrementing 'dp += len' speed went from 70 to 45  (167 to 116 for double)
	 */
        i = cnt;
	for (i=0, dp=data; i<cnt; i++, dp += len) {
	  if (Qmemcpy) {
	    /* 40 for double */
	    byteswap_doubles((double *)dp);
	  } else
	    /* 116 for double */
	    bswap(dp,len,1);
	}
      }
    } while (repeat-- > 0);
    nrepeat++;
    t2 = 60 * cputime();
    if (t1==t0) warning("Inaccurate measurement for reading");
    else rspeed=cnt*len/(t1-t0)/1048576.0*nrepeat;
    if (t1==t2) warning("Inaccurate measurement for writing");
    else wspeed=cnt*len/(t2-t1)/1048576.0*nrepeat;
    dprintf(1,"Read %d: %f s; %d x swap %d * %d bytes: %f s; %g %g Mswap\n",
       		cnt*len, t1-t0, nrepeat, cnt, len, t2-t1, 
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

