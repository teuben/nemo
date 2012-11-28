/*
 * FILESECRET.C: binary file structure package.
 *
 *  Warning: ZENO's version is labeled version 3, our version 3
 *           has nothing to do with their version 3.  There are some
 *           minor differences in the implementation that result in
 *           a difference in capabilities, e.g. ZENO can handle
 *           multiple same-named items within a set, whereas NEMO
 *           cannot
 *
 *  Input is done through:      Output through:
 *      fread()                     fwrite()
 *      getxstr() -> getc()         putxstr() ->  putc()    (see extstring.c)
 *      saferead -> fread()         fwrite()
 *
 * V 1.0: Joshua Barnes 4/87	basic I/O operators implemented,
 * V 1.X: Lyman Hurd 8/87	adde f-d coercion, deferred input,
 * V 2.0: Joshua Barnes 4/88	new types, operators, external format.
 * V 2.1: Peter Teuben:	4/90	bug in stacklength in filesecret.h removed
 * V 2.2: Peter Teuben: 9/90    random access I/O - swap
 * V 2.4: Peter Teuben: 3/91    type recasting, essentially for csf.c
 * V 2.5: 13-feb-92   PJT	added bool qsf() to query if a stru.file?
 *        10-nov-92   PJT	pretty bad bug in f2d/d2f conversion routines
 * V 2.6a 22-feb-94   PJT       ansi
 *      b 24-jun-94   PJT       implemented float_to_double conversiom (csf)
 *      c 20-nov-94   PJT       added <unistd.h>
 *	  12-apr-05   pjt       use NOPROTO instead of PROTO
 *      d 16-feb-97   pjt       extra protection to malloc(<0)
 *      e 13-jan-99   pjt       defer free in random I/O
 *      f  6-apr-01   pjt       malloc->calloc
 *      g  7-jun-01   pjt       aptr -> ap (why did the compiler not warn?)
 * V 2.7  20-jun-01   pjt       no more NOPROTO, fixed protos for gcc3
 *      a  8-oct-01   pjt       flush buffer when put_tes() at top level
 * V 3.0  23-may-02   pjt    Support for >2GB (large file size)
 *        17-mar-03   pjt    fix serious memory usage bug for f2d conversion
 *        18-jun-03   pjt    <docs>
 * V 3.1  15-mar-05   pjt    C++ compilable
 * V 3.2   2-jun-05   pjt    blocked (sequential) I/O
 * V 3.3  25-may-07   pjt    handle > 2GB objects in memory (Pierre Fortin <pierre.fortin@oamp.fr>)
 * V 3.4  12-dec-09   pjt    support the new halfp type for I/O (see also csf)
 *        27-Sep-10   jcl    MINGW32/WINDOWS support
 *
 *  Although the SWAP test is done on input for every item - for deferred
 *  input it may fail if in the mean time another file was read which was
 *  not in swap mode => Perhaps better to HardCode (#define)???
 *
 *  array of strings ??
 *                                  
 */

#include <stdinc.h>
#include <unistd.h>
#include <strlib.h>
#include <filestruct.h>
#include <extstring.h>
#include "filesecret.h"
#include <stdarg.h>

extern int convert_d2f(int, double *, float  *);
extern int convert_f2d(int, float  *, double *);
extern int convert_h2f(int, halfp  *, float  *);
extern int convert_h2d(int, halfp  *, double *);
extern int convert_f2h(int, float  *, halfp  *);
extern int convert_d2h(int, double *, halfp  *);
#ifdef __MINGW32__
#define fseeko fseek
#define ftello ftell
#endif

/*
 * COPY_ITEM: recursively copy item from input to output.
 * An example of recursive file traversal and memory etiquette.
 */

void copy_item(stream ostr, stream istr, string tag)
{
    string type, *tags, *tp;
    int *dims;
    size_t dlen;
    byte *buf;

    if (! get_tag_ok(istr, tag))		/* prevent obvious errors   */
	error("copy_item: tag %s not found", tag);
    type = get_type(istr, tag);			/* find out type of data    */
    if (! streq(type, SetType)) {		/* a basic type or array?   */
	dims = get_dims(istr, tag);		/*   find out about shape   */
	dlen = get_dlen(istr, tag);		/*   and length in bytes    */
	if(dlen<0) error("copy_item: %s with dlen=%d",tag,dlen);    /* yuck */
	buf = (byte *) calloc(dlen,1);		/*   get space for a buffer */
	if (buf == NULL)			/*   and check for error    */
	    error("copy_item: item %s: not enuf memory", tag);
	get_data_sub(istr, tag, type, buf, dims, FALSE);
						/*   read data from input   */
	put_data_sub(ostr, tag, type, buf, dims, FALSE);
						/*   and write it to output */
	if (dims != NULL)			/*   free dimension list    */
	    free(dims);
	free(buf);				/*   free temporary buffer  */
    } else {					/* a set of other items?    */
	get_set(istr, tag);			/*   access set's contents  */
	put_set(ostr, tag);			/*   output set token       */
	tags = list_tags(istr);			/*   get list of set tags   */
	for (tp = tags; *tp != NULL; tp++)	/*   loop tp over them      */
	    copy_item(ostr, istr, *tp);		/*     copy each in turn    */
	get_tes(istr, tag);			/*   close access to set    */
	put_tes(ostr, tag);			/*   output termial symbol  */
	for (tp = tags; *tp != NULL; tp++)	/*   loop over tags again   */
	    free(*tp);				/*     free each one up     */
	free(tags);				/*   free tag vector itself */
    }
    free(type);					/* free up type string      */
}
/*
 * COPY_ITEM_CVT: recursively copy item from input to output.
 *	With the option to convert data type
 */

void copy_item_cvt(stream ostr, stream istr, string tag, string *cvt)
{
    string type, *tags, *tp;
    int *dims, cvtlen;
    size_t dlen;
    char *cp;
    byte *bufin, *bufout=NULL;
    itemptr ipt;

    if (! get_tag_ok(istr, tag))		/* prevent obvious errors   */
	error("copy_item_cvt: tag %s not found", tag);
    type = get_type(istr, tag);			/* find out type of data    */
    cvtlen = xstrlen(cvt,sizeof(string)) - 1;   /* # conversions            */
    if (! streq(type, SetType)) {		/* a basic type or array?   */
	dims = get_dims(istr, tag);		/*   find out about shape   */
	dlen = get_dlen(istr, tag);		/*   and length in bytes    */
	if(dlen<0) error("copy_item_cvt: %s with dlen=%d",tag,dlen); 
	bufin = (byte *) calloc(dlen,1);	/*   get space for a buffer */
	if (bufin == NULL)			/*   and check for error    */
	    error("copy_item_cvt: item %s: not enuf memory", tag);
	get_data_sub(istr, tag, type, bufin, dims, FALSE);
						/*   read data from input   */
        cp = findtype(cvt,type);                /* check if to convert      */
        if (cp==NULL)
	    put_data_sub(ostr, tag, type, bufin,  dims, FALSE); /* copy now */
        else if (streq(type,DoubleType)) {
            if (streq(cp,"d2f")) {		/* convert double to float */
                dprintf(1,"Converting %s in %s\n",cp,tag);
                ipt = makeitem(FloatType,tag,NULL,dims);    /* silly */
                convert_d2f(eltcnt(ipt,0),(double*)bufin,(float*)bufin);
	        put_data_sub(ostr, tag, FloatType, bufin,  dims, FALSE); 
                freeitem(ipt,FALSE);
	    } else if (streq(cp,"d2h")) {
                dprintf(1,"Converting %s in %s\n",cp,tag);
                ipt = makeitem(HalfpType,tag,NULL,dims);    /* silly */
                convert_d2h(eltcnt(ipt,0),(double*)bufin,(halfp*)bufin);
	        put_data_sub(ostr, tag, HalfpType, bufin,  dims, FALSE); 
                freeitem(ipt,FALSE);
            } else {
            	warning("Cannot convert %s yet in %s",cp,tag);
	        put_data_sub(ostr, tag, type, bufin,  dims, FALSE); 
	    }
	} else if (streq(type,FloatType)) {
            if (streq(cp,"f2d")) {		/* convert float to double */
                dprintf(1,"Converting %s in %s\n",cp,tag);
                ipt = makeitem(DoubleType,tag,NULL,dims);    /* silly */
                bufout = (byte *) allocate(datlen(ipt,0));
                if (bufout == NULL)
               	    error("copy_item_cvt: item %s: (f2d) not enuf memory", tag);
                convert_f2d(eltcnt(ipt,0),(float*)bufin,(double*)bufout);
	        put_data_sub(ostr, tag, DoubleType, bufout,  dims, FALSE); 
                freeitem(ipt,0);
	    } else if (streq(cp,"f2h")) {
                dprintf(1,"Converting %s in %s\n",cp,tag);
                ipt = makeitem(HalfpType,tag,NULL,dims);    /* silly */
                bufout = (byte *) allocate(datlen(ipt,0));
                if (bufout == NULL)
               	    error("copy_item_cvt: item %s: (f2h) not enuf memory", tag);
                convert_f2h(eltcnt(ipt,0),(float*)bufin,(halfp*)bufout);
	        put_data_sub(ostr, tag, HalfpType, bufout,  dims, FALSE); 
                freeitem(ipt,0);
            } else {
            	warning("Cannot convert %s yet in %s",cp,tag);
	        put_data_sub(ostr, tag, type, bufin,  dims, FALSE); 

	    }
	} else if (streq(type,HalfpType)) {
            if (streq(cp,"h2d")) {		/* convert halfp to double */
                dprintf(1,"Converting %s in %s\n",cp,tag);
                ipt = makeitem(DoubleType,tag,NULL,dims);    /* silly */
                bufout = (byte *) allocate(datlen(ipt,0));
                if (bufout == NULL)
               	    error("copy_item_cvt: item %s: (h2d) not enuf memory", tag);
                convert_h2d(eltcnt(ipt,0),(halfp*)bufin,(double*)bufout);
	        put_data_sub(ostr, tag, DoubleType, bufout,  dims, FALSE); 
                freeitem(ipt,0);
	    } else if (streq(cp,"h2f")) {
                dprintf(1,"Converting %s in %s\n",cp,tag);
                ipt = makeitem(FloatType,tag,NULL,dims);    /* silly */
                bufout = (byte *) allocate(datlen(ipt,0));
                if (bufout == NULL)
               	    error("copy_item_cvt: item %s: (h2f) not enuf memory", tag);
                convert_h2f(eltcnt(ipt,0),(halfp*)bufin,(float*)bufout);
	        put_data_sub(ostr, tag, FloatType, bufout,  dims, FALSE); 
                freeitem(ipt,0);
            } else {
            	warning("Cannot convert %s yet in %s",cp,tag);
	        put_data_sub(ostr, tag, type, bufin,  dims, FALSE); 
	    }
	} else if (streq(type,IntType)) {
            warning("Cannot convert %s yet in %s",cp,tag);
	    put_data_sub(ostr, tag, type, bufin,  dims, FALSE); /* copy now */
	} else if (streq(type,ShortType)) {
            warning("Cannot convert %s yet in %s",cp,tag);
	    put_data_sub(ostr, tag, type, bufin,  dims, FALSE); /* copy now */
	} else {
	    if (cvtlen>0) warning("Cannot convert type %c in %s",type,tag);
	    put_data_sub(ostr, tag, type, bufin,  dims, FALSE); /* copy now */
        }
						/*   and write it to output */
	if (dims != NULL)			/*   free dimension list    */
	    free(dims);
	free(bufin);		                /*   free temporary buffer  */
        if (bufout) free(bufout);               /* if used: free this too   */
    } else {					/* a set of other items?    */
	get_set(istr, tag);			/*   access set's contents  */
	put_set(ostr, tag);			/*   output set token       */
	tags = list_tags(istr);			/*   get list of set tags   */
	for (tp = tags; *tp != NULL; tp++)	/*   loop tp over them      */
	    copy_item_cvt(ostr, istr, *tp, cvt);/*     copy each in turn    */
	get_tes(istr, tag);			/*   close access to set    */
	put_tes(ostr, tag);			/*   output termial symbol  */
	for (tp = tags; *tp != NULL; tp++)	/*   loop over tags again   */
	    free(*tp);				/*     free each one up     */
	free(tags);				/*   free tag vector itself */
    }
    free(type);					/* free up type string      */
}

/************************************************************************/
/*                         USER OUTPUT FUNCTIONS                        */
/************************************************************************/

/*
 * PUT_SET: begin named set in output.
 */

void put_set(stream str, string tag)
{
    strstkptr sspt;
    itemptr ipt;

    sspt = findstream(str);			/* get stream-stack struct  */
    ipt = makeitem(SetType, tag, NULL, NULL);	/* make item to hold tag    */
    ss_push(sspt, ipt);				/* and stack for put_tes    */
    put_data(str, tag, SetType, NULL, 0);	/* output external token    */
}  

/*
 * PUT_TES: end set in output.
 */

void put_tes(stream str, string tag)
{
    strstkptr sspt;
    itemptr ipt;

    sspt = findstream(str);			/* get stream-stack struct  */
    if (sspt->ss_stp < 0)			/* check for underflow      */
	error("put_tes: stack underflow");
    ipt = sspt->ss_stk[sspt->ss_stp];		/* get top item on stack    */
    if (tag != NULL &&				/* tag is specified	    */
	  ! streq(ItemTag(ipt), tag))		/*   but does not match?    */
	error("put_tes: set = %s tes = %s", ItemTag(ipt), tag);
    sspt->ss_stk[sspt->ss_stp] = NULL;		/* erase pointer to item    */
    freeitem(ipt, FALSE);			/* and reclaim storage      */
    ss_pop(sspt);				/* flush stacked item       */
    put_data(str, NULL, TesType, NULL, 0);	/* output external token    */
    if (sspt->ss_stp == -1) {                   /* if at top level          */
      dprintf(1,"put_tes(%s) flushing\n",tag);  /* removed '\n' 27/06/08 WD */
      fflush(str);                              /* flush buffer for Walter  */ 
    }
}

/*
 * PUT_STRING: write string to a structured file.
 */

void put_string(stream str, string tag, string dat)
{
    put_data(str, tag, CharType, dat, xstrlen(dat, 1), 0);
}

/*
 * PUT_DATA: write data object to a structured file.
 * Synopsis: put_data(str, tag, typ, dat, dimN, ..., dim1, 0)
 */
void put_data(stream str, string tag, string typ, void *dat, int dim1, ...)
{
    va_list ap;
    int dim[MaxVecDim], n = 0;

    dim[0] = dim1;

    va_start(ap, dim1);				/* access argument list     */
    while (dim[n++] > 0) {			/* loop reading dimensions  */
	if (n >= MaxVecDim)			/*   no room for any more?  */
	    error("put_data: too many dims; item %s", tag);
	dim[n] = va_arg(ap, int);		/*   else get next argument */
    }
    va_end(ap);
    						/* call next level routine  */
    put_data_sub(str, tag, typ, dat, (dim[0] != 0 ? dim : NULL), FALSE);
}
/*
 * PUT_DATA_SUB: worker for above manager.
 */
void put_data_sub(
    stream str, 	/* stream to write data to */
    string tag, 	/* tag for output item */
    string typ,     	/* data type for output */
    void *dat,	 	/* place to store data */
    int *dim,		/* vector of dimensions */
    bool con)		/* coercion flag (not used) */
{
    itemptr ipt;

    ipt = makeitem(typ, tag, dat, dim);		/* make item wo/ copying    */
    if (! putitem(str, ipt)) 			/* output external rep.     */
	error("put_data_sub: putitem failed");
    freeitem(ipt, FALSE);			/* and reclaim storage      */
}

/************************************************************************/
/*                         USER OUTPUT FUNCTIONS (RANDOM)               */
/************************************************************************/
#if defined(RANDOM)
/*
 * PUT_DATA_SET: open an item for random access
 * Synopsis: put_data_set(str, tag, typ, dimN, ..., dim1, 0)
 */
void put_data_set(stream str, string tag, string typ, int dim1, ...)
{
    va_list ap;
    int dim[MaxVecDim], *buf, n = 0;
    itemptr ipt;
    strstkptr sspt;

    dim[0] = dim1;
    va_start(ap, dim1);				/* access argument list     */
    while (dim[n++] > 0) {			/* loop reading dimensions  */
	if (n >= MaxVecDim)			/*   no room for any more?  */
	    error("put_data_set: too many dims; item %s", tag);
	dim[n] = va_arg(ap, int);		/*   else get next argument */
    } while (dim[n++] != 0);			/* until a zero comes up    */
    va_end(ap);

    sspt = findstream(str);
    if (sspt->ss_ran)
        error("put_data_set: %s: can currently handle one random access item",tag);
    buf = (int *) copxstr(dim,sizeof(int));
    ipt = makeitem(typ,tag,NULL,buf);            /* make item but no copy */
    sspt->ss_ran = ipt;
    puthdr(str,ipt);                           /* write the header right now */

    ItemPos(ipt) = ftello(str);                      /* begin of random data */
    ItemOff(ipt) = 0;                                   /* offset where I/O */
    sspt->ss_pos = ftello(str) + datlen(ipt,0);         /* end of random data */
}

/*
 * PUT_DATA_TES: close an item for random access
 * Synopsis:  put_data_tes(str, tag)
 */

void put_data_tes(stream str, string tag)
{
    itemptr ipt;
    strstkptr sspt;

    sspt = findstream(str);
    ipt = sspt->ss_ran;
    if (ipt == NULL) error("put_data_tes: item %s is not random",tag);
    ipt = sspt->ss_ran;
    if (!streq(tag,ItemTag(ipt))) error("put_data_tes: invalid tag name %s",tag);
    fseeko(str,sspt->ss_pos,0);      /* go where we left off */
    sspt->ss_pos = 0L;              /* mark file i/o sequential again */
    sspt->ss_ran = NULL;            /* reset pointer to the random item */
    free( ItemDim(ipt));
    freeitem(ipt,FALSE);
}
/*
 * PUT_DATA_RAN: random acces output of data
 * Synopsis:   put_data_ran(str, tag, dat, offset, length)
 */

void put_data_ran(
    stream str,
    string tag,
    void *dat,
    int offset,
    int length)
{
    itemptr ipt;
    strstkptr sspt;

    sspt = findstream(str);
    ipt = sspt->ss_ran;
    if (ipt==NULL) error("put_data_ran: tag %s no random item",tag);
    if (!streq(tag,ItemTag(ipt))) error("put_data_ran: invalid tag name %s",tag);
    offset *= ItemLen(ipt);     /* as input offset and length were  */
    length *= ItemLen(ipt);     /* in units of itemlen !!! */
    if (offset+length > datlen(ipt,0))
        error("put_data_ran: tag %s cannot write beyond allocated boundary",tag);
    fseeko(str,offset + ItemPos(ipt),0);
    if (length != fwrite((char *)dat,sizeof(byte),length,str))
        error("put_data_ran: error writing tag %s",tag);
}

/*
 * PUT_DATA_BLOCKED: blocked sequential acces output of data
 * Synopsis:   put_data_blocked(str, tag, dat, length)
 */

void put_data_blocked(
    stream str,
    string tag,
    void *dat,
    int length)
{
    itemptr ipt;
    strstkptr sspt;
    int offset;

    sspt = findstream(str);
    ipt = sspt->ss_ran;
    if (ipt==NULL) error("put_data_blocked: tag %s no random item",tag);
    if (!streq(tag,ItemTag(ipt))) error("put_data_blocked: invalid tag name %s",tag);
    offset = ItemOff(ipt);
    length *= ItemLen(ipt);     /* in units of itemlen !!! */
    if (offset+length > datlen(ipt,0))
        error("put_data_blocked: tag %s cannot write beyond allocated boundary",tag);
    // no fseek() needed in blocked() !!!!
    // fseeko(str,offset + ItemPos(ipt),0);
    if (length != fwrite((char *)dat,sizeof(byte),length,str))
        error("put_data_blocked: error writing tag %s",tag);
    ItemOff(ipt) += length;
}

#else

    /* Un-implemented parts of random access write */

void put_data_set(va_alist)
va_dcl
{
    va_list aptr;
    string tag;

    va_start(aptr);				/* access argument list     */
    (void) va_arg(aptr, stream);		/* collect mandatory args   */
    tag = va_arg(aptr, string);
    warning("put_data_set: skipped %s no random access",tag);
}

void put_data_tes()
{
}

void put_data_ran()
{
}
#endif

/************************************************************************/
/*                          USER INPUT FUNCTIONS                        */
/************************************************************************/

/*
 * GET_SET: access named set in input.
 */

void get_set(stream str, string tag)
{
    strstkptr sspt;
    itemptr sptr;
  
    sspt = findstream(str);			/* get associated info      */
    sptr = scantag(sspt, tag);			/* scan for that tag        */
    if (sptr == NULL)				/* check for end of file    */
	error("get_set: at EOF");
    if (! streq(ItemTyp(sptr), SetType))	/* be sure its a set        */
	error("get_set: %s not a set", tag);
    ss_push(sspt, sptr);
}

/*
 * GET_TES: end access of named set.
 */

void get_tes(stream str, string tag)
{
    itemptr ipt;
    strstkptr sspt;

    sspt = findstream(str);			/* get associated info      */
    if (sspt->ss_stp < 0)			/* catch stack underflow    */
	error("get_tes: stream stack underflow");
    ipt = sspt->ss_stk[sspt->ss_stp];		/* get item on top of stack */
    if (tag != NULL &&				/* is tag specified, but    */
	  ! streq(ItemTag(ipt), tag))		/* does not match?          */
	error("get_tes: set = %s tes = %s", ItemTag(ipt), tag);
    ss_pop(sspt);				/* remove item from stack   */
    if (sspt->ss_stp == -1) {			/* back to top level?	    */
	freeitem(sspt->ss_stk[0], TRUE);	/*   then free input set    */
	sspt->ss_stk[0] = NULL;			/*   and flush pending item */
    }
}

/*
 * GET_DATA: read data object from a structured file.
 * Synopsis: get_data(str, tag, typ, dat, dimN, ..., dim1, 0)
 */
void get_data(stream str, string tag, string typ, void *dat, int dim1, ...)
{
    va_list ap;
    int dim[MaxVecDim], n = 0;

    dim[0] = dim1;
    va_start(ap, dim1);				/* access argument list     */
    while (dim[n++] > 0) {			/* loop reading dimensions  */
	if (n >= MaxVecDim)			/*   no room for any more?  */
	    error("get_data: item %s: too many dims", tag);
	dim[n] = va_arg(ap, int);		/*   else get next argument */
    } 
    va_end(ap);
    						/* call next level routine  */
    get_data_sub(str, tag, typ, dat, (dim[0] != 0 ? dim : NULL), FALSE);
}
/*
 * GET_DATA_COERCED : read data object from a structured file,
 * automatically performing float <--> double conversions.
 * Synopsis: get_data_coerced(str, tag, typ, dat, dimN, ..., dim1, 0)
 */
void get_data_coerced(stream str,string tag,string typ,void *dat,int dim1, ...)
{
    va_list ap;
    int dim[MaxVecDim], n = 0;

    dim[0] = dim1;
    va_start(ap,dim1);				/* access argument list     */
    while (dim[n++] > 0) {			/* loop reading dimensions  */
	if (n >= MaxVecDim)			/*   no room for any more?  */
	    error("get_data_coerced: item %s: too many dims", tag);
	dim[n] = va_arg(ap, int);		/*   else get next argument */
    } 
    va_end(ap);
    						/* call next level routine  */
    get_data_sub(str, tag, typ, dat, (dim[0] != 0 ? dim : NULL), TRUE);
}
	
/*
 * GET_DATA_SUB: worker for above managers.
 */
void get_data_sub(
    stream str,             	/* stream to read data from */
    string tag,             	/* expected item tag */
    string typ,               	/* expected data type */
    void *dat,              	/* place to store data */
    int *dim,			/* array of dimensions */
    bool con)			/* coercion flag */
{
    strstkptr sspt;
    itemptr ipt;
    copyproc cop;

    sspt = findstream(str);			/* access assoc. info	    */
    ipt = scantag(sspt, tag);			/* scan input for tag	    */
    if (ipt == NULL)				/* check input succeeded    */
	error("get_data: at EOF");
    if (! con) {				/* rigid type checking?     */
	if(! streq(typ, ItemTyp(ipt)))		/*   and types dont match?  */
	    error("get_data_sub: item %s: types %s, %s don't match",
		  tag, ItemTyp(ipt), typ);
	cop = copydata;				/*   get default routine    */ /*C++*/
    } else {
	cop = copyfun(ItemTyp(ipt), typ);	/*   get specialist routine */
	if (cop == NULL)			/*   but be sure one exists */
	    error("get_data_sub: item %s: types %d, %d don't convert",
		  tag, ItemTyp(ipt), typ);
    }
    if (dim != NULL && ItemDim(ipt) != NULL &&	/* check layout of data     */
	  ! xstreq(dim, ItemDim(ipt), sizeof(int)))
	error("get_data_sub: item %s: dimensions don't match", tag);
    else if (dim == NULL && ItemDim(ipt) != NULL)
	error("get_data_sub: item %s: can't copy plural to scalar", tag);
    else if (dim != NULL && ItemDim(ipt) == NULL)
	error("get_data_sub: item %s: can't copy scalar to plural", tag);
    (cop)(dat, 0, eltcnt(ipt,0), ipt, str);    	/* copy data from input     */ /*C++*/
    if (sspt->ss_stp == -1)			/* was input at top level?  */
	freeitem(ipt, TRUE);			/*   yes, free saved item   */
}

/************************************************************************/
/*                          USER INPUT FUNCTIONS (RANDOM)               */
/************************************************************************/
#if defined(RANDOM)
/*
 * GET_DATA_SET: open an item for random access
 * Synopsis: get_data_set(str, tag, typ, dimN, ..., dim1, 0)
 */
void get_data_set(stream str, string tag, string typ, int dim1, ...)
{
    va_list ap;
    int dim[MaxVecDim], n = 0;
    itemptr ipt;
    strstkptr sspt;

    dim[0] = dim1;

    va_start(ap, dim1);				/* access argument list     */
    while (dim[n++] > 0) {			/* loop reading dimensions  */
	if (n >= MaxVecDim)			/*   no room for any more?  */
	    error("put_data_set: too many dims; item %s", tag);
	dim[n] = va_arg(ap, int);		/*   else get next argument */
    } 
    va_end(ap);

    sspt = findstream(str);
    if (sspt->ss_ran)
        error("put_data_set: %s: can only handle one random access item",tag);
    ipt = scantag(sspt,tag);	/* try and find the data */
    if (ipt==NULL) error("get_data_set: Bad EOF");

    sspt->ss_pos = ItemPos(ipt) + datlen(ipt,0);         /* end of random data */
    sspt->ss_ran = ipt;
}

/*
 * GET_DATA_TES: close an item for random access
 * Synopsis:  get_data_tes(str, tag)
 */

void get_data_tes(stream str, string tag)
{
    itemptr ipt;
    strstkptr sspt;

    sspt = findstream(str);
    ipt = sspt->ss_ran;
    if (ipt == NULL) error("get_data_tes: item %s is not random",tag);
    ipt = sspt->ss_ran;
    if (!streq(tag,ItemTag(ipt))) error("get_data_tes: invalid tag name %s",tag);

    sspt->ss_pos = 0L;              /* mark file i/o sequential again */
    sspt->ss_ran = NULL;            /* reset pointer to the random item */
#if 0
	/* should delay this */
    free( ItemDim(ipt));
    freeitem(ipt,FALSE);
#endif
}

void get_data_ran(
    stream str,
    string tag,
    void *dat,
    int offset,
    int length
) {
    itemptr ipt;
    strstkptr sspt;

    sspt = findstream(str);
    ipt = sspt->ss_ran;
    if (ipt==NULL)
        error("get_data_ran: tag %s is not in random access mode",tag);
    copydata(dat,offset,length,ipt,str);
}

void get_data_blocked(
    stream str,
    string tag,
    void *dat,
    int length
) {
    itemptr ipt;
    strstkptr sspt;
    int offset;

    sspt = findstream(str);
    ipt = sspt->ss_ran;
    offset = ItemOff(ipt);
    if (ipt==NULL)
        error("get_data_blocked: tag %s is not in blocked access mode",tag);
    copydata(dat,offset,length,ipt,str);
    ItemOff(ipt) = offset+length;
}

#endif


/*
 * GET_STRING: read a string from a structured file.
 *             returns pointer to freshly allocated space that can be free'd
 */

string get_string(
    stream str,			/* stream to read string from */
    string tag			/* tag to compare with item */
) {
    strstkptr sspt;
    itemptr ipt;
    int *dp;
    size_t dlen;
    char *dat;

    sspt = findstream(str);			/* access assoc. info	    */
    ipt = scantag(sspt, tag);			/* scan input for tag	    */
    if (ipt == NULL)				/* check input succeeded    */
	error("get_string: at EOF");
    dp = ItemDim(ipt);				/* get list of dimensions   */
    if (! streq(ItemTyp(ipt), CharType) ||	/* check type of item       */
	  dp == NULL || *dp++ == 0 || *dp != 0)	/* and shape of data        */
	error("get_string: item %s: not plural char", tag);
    dlen = datlen(ipt,0);
    if(dlen<0) error("get_string: %s with dlen=%d",tag,dlen);       /* yuck */
    dat = (char *) calloc(dlen,1);        	/* alloc memory for string  */
    if (dat == NULL)				/* did alloc fail?          */
	error("get_string: item %s: not enuf memory", tag);
    copydata(dat, 0, dlen, ipt, str);		/* copy string from input   */
    if (sspt->ss_stp == -1)			/* item read at top level?  */
	freeitem(ipt, TRUE);			/*   yes, so free it up     */
    return (dat);				/* return string as value   */
}

/*
 * GET_TAG_OK: determine if subsequent get_data(), get_set(),
 * operation will find specified tag.  Returns FALSE at EOF.
 */

bool get_tag_ok(
    stream str,			/* input stream obtained from stropen */
    string tag			/* tag of item being looked for */
) {
    strstkptr sspt;
    itemptr ipt;

    sspt = findstream(str);			/* lookup associated entry  */
    if (sspt->ss_stp == -1) {			/* input from top level?    */
	ipt = nextitem(sspt);			/*   look at next item      */
	return (ipt != NULL && streq(tag, ItemTag(ipt)));
						/*   test existance and tag */
    } else					/* input within a set?      */
	return (finditem(sspt, tag) != NULL);	/*   test success of scan   */
}

/*
 * LIST_TAGS: list all tags for which get_tag_ok is true.
 */

string *list_tags(
stream str
) {
    strstkptr sspt;
    itemptr ipt, *spt;
    string tags[MaxSetLen], *tpt = &tags[0];

    sspt = findstream(str);			/* lookup associated entry  */
    if (sspt->ss_stp == -1) {			/* input from top level?    */
	ipt = nextitem(sspt);			/*   get next item read in  */
	if (ipt == NULL)			/*   nothing left in input? */
	    return (NULL);			/*     then return nothing  */
	*tpt++ = (string) copxstr(ItemTag(ipt), sizeof(char));
						/*   make copy of item tag  */
    } else {					/* input within a set?      */
	ipt = sspt->ss_stk[sspt->ss_stp];	/*   get current item set   */
	spt = (itemptr *) ItemDat(ipt);		/*   get string of items    */
	while (*spt != NULL)			/*   loop over items        */
	    *tpt++ = (string) copxstr(ItemTag(*spt++), sizeof(char));
						/*   make copy of item tag  */
    }
    *tpt = NULL;				/* terminate string of tags */
    return ((string *) copxstr(tags, sizeof(string)));
    						/* return copy of copies    */
}	

/*
 * GET_TYPE, GET_DIMS, GET_DLEN: enquire about type, dimensions, and data
 * length of specified item.
 */

string get_type(
    stream str,			/* input stream obtained from stropen */
    string tag			/* tag of item being looked for */
) {
    strstkptr sspt;
    itemptr ipt;

    sspt = findstream(str);			/* lookup associated entry  */
    ipt = scantag(sspt, tag);			/* and scan for named item  */
    if (ipt == NULL)				/* check scan succeeded     */
	error("get_type: at EOF");
    if (sspt->ss_stp == -1)			/* was input at top level?  */
	sspt->ss_stk[0] = ipt;			/*   put back for next time */
    return ((string) copxstr(ItemTyp(ipt), sizeof(char)));
						/* return copy of type      */
}

int *get_dims(
    stream str,			/* input stream obtained from stropen */
    string tag			/* tag of item being looked for */
) {
    strstkptr sspt;
    itemptr ipt;

    sspt = findstream(str);			/* lookup associated entry  */
    ipt = scantag(sspt, tag);			/* and scan for named item  */
    if (ipt == NULL)				/* check scan succeeded     */
	error("get_dims: at EOF");
    if (sspt->ss_stp == -1)			/* was input at top level?  */
	sspt->ss_stk[0] = ipt;			/*   put back for next time */
    if (ItemDim(ipt) != NULL)			/* any dimensions to item?  */
	return ((int *) copxstr(ItemDim(ipt), sizeof(int)));
						/*   return copy of dims    */
    else
	return NULL;
}

size_t get_dlen(
    stream str,			/* input stream obtained from stropen */
    string tag			/* tag of item being looked for */
) {
    strstkptr sspt;
    itemptr ipt;

    sspt = findstream(str);			/* lookup associated entry  */
    ipt = scantag(sspt, tag);			/* and scan for named item  */
    if (ipt == NULL)				/* check scan succeeded     */
	error("get_dlen: at EOF");
    if (sspt->ss_stp == -1)			/* was input at top level?  */
	sspt->ss_stk[0] = ipt;			/*   put back for next time */
    return datlen(ipt, 0);			/* return count of bytes    */
}

/*
 * SKIP_ITEM: flush pending item in top level input, nop within set.
 * Returns FALSE at EOF.
 */

bool skip_item(
    stream str			/* input stream obtained from stropen */
) {
    strstkptr sspt;
    itemptr ipt;

    sspt = findstream(str);			/* lookup associated entry  */
    if (sspt->ss_stp == -1) {			/* input from top level?    */
	ipt = nextitem(sspt);			/*   get item read next     */
	if (ipt == NULL)			/*   nothing left in input? */
	    return FALSE;			/*     then nothing to skip */
	freeitem(ipt, TRUE);			/*   reclaim storage space  */
	sspt->ss_stk[0] = NULL;			/*   and flush that item    */
	return TRUE;				/*   handled an item        */
    } else {
	printf("skip_item: within set");
	return TRUE;				/*    a lie                 */
    }
}

/************************************************************************/
/*                                OUTPUT                                */
/************************************************************************/

/*
 * WRITEITEM: write an item or set of items.     [NOT USED]
 */

local bool writeitem(stream str, itemptr ipt)
{
    itemptr *setp, tesp;

    if (! streq(ItemTyp(ipt), SetType))		/* just a simple item?      */
	return (putitem(str, ipt));		/*   then write it out      */
    else {					/* recursively handle set   */
	if (! putitem(str, ipt))		/*   write out set item     */
	    return (FALSE);
	setp = (itemptr*) ItemDat(ipt);		/*   point to item string   */
	while (*setp != NULL)			/*   loop over the items    */
	    writeitem(str, *setp++);		/*     write each one out   */
	tesp = makeitem(TesType, NULL, NULL, NULL);
						/*   create closing item    */
	if (! putitem(str, tesp))		/*   write out tes item     */
	    return (FALSE);
	freeitem(tesp, FALSE);			/*   and reclaim memory     */
	return (TRUE);				/*   indicate success       */
    }
}

/*
 * PUTITEM: write item to stream; returns indication of success/failure.
 */

local bool putitem(stream str, itemptr ipt)
{
    if (! puthdr(str, ipt))                     /* write item header        */
        return (FALSE);
    if (! streq(ItemTyp(ipt), SetType) && ! streq(ItemTyp(ipt), TesType))
						/* an ordinary data item?   */
	if (! putdat(str, ipt))                 /*   write item data        */
	    return (FALSE);
    return (TRUE);                              /* indicate success         */
}

/*
 * PUTHDR: write item header to stream.
 */

local bool puthdr(stream str, itemptr ipt)
{
    short num;

    
    num = (ItemDim(ipt) == NULL) ? SingMagic : PlurMagic;
    						/* determine magic number   */
    if (fwrite((char *)&num, sizeof(short), 1, str) != 1)
	return (FALSE);				/* return FALSE on failure  */
    if (! putxstr(str, ItemTyp(ipt), sizeof(char)))
	return (FALSE);
    if (ItemTag(ipt) != NULL) {                 /* is item tagged?          */
        if (xstrlen(ItemTag(ipt), sizeof(char)) > MaxTagLen)
            error("puthdr: tag too long");
        if (! putxstr(str, ItemTag(ipt), sizeof(char)))
	    return (FALSE);                     /*   write item tag         */
    }
    if (ItemDim(ipt) != NULL) {                 /* a vectorized item?       */
        if (xstrlen(ItemDim(ipt), sizeof(int)) > MaxVecDim)
            error("puthdr: too many dimensions");
        if (! putxstr(str, ItemDim(ipt), sizeof(int)))
	    return (FALSE);                     /*   write vect dims        */
    }
    return(TRUE);                               /* indicate success         */
}

/*
 * PUTDAT: write data of a item.
 */

local bool putdat(stream str, itemptr ipt)
{
    size_t len;

    if (ItemDat(ipt) == NULL)			/* no data to write?        */
	error("putdat: item %s has no data", ItemTag(ipt));
    len = datlen(ipt, 0);			/* count bytes to output  */
    return (fwrite((char*)ItemDat(ipt), sizeof(byte), len, str) == len);
						/* write data to stream   */
}

/************************************************************************/
/*                                 INPUT                                */
/************************************************************************/

/*
 * SCANTAG: look for tagged item in input.  Returns NULL on EOF,
 * calls error() if requested item not found.
 */

local itemptr scantag(strstkptr sspt, string tag)
{
    itemptr ipt;

    if (sspt->ss_stp == -1) {			/* input at top-level?    */
	ipt = nextitem(sspt);			/*   check out next item  */
	sspt->ss_stk[0] = NULL;			/*   take charge of it    */
	if (ipt != NULL && ! streq(tag, ItemTag(ipt)))
	    					/*   check tags match	  */
	    error("scantag: got %s instead of %s",
		  ItemTag(ipt), tag);
    } else {					/* input within a set?	  */
	ipt = finditem(sspt, tag);		/*   scan for the tag     */
	if (ipt == NULL)			/*   was none found?	  */
	    error("scantag: item %s not found in set %s",
		  tag, ItemTag(sspt->ss_stk[sspt->ss_stp]));
    }
    return (ipt);				/* return item found	  */
}

local itemptr nextitem(strstkptr sspt)
{
    itemptr ipt;

    if (sspt->ss_stk[0] != NULL)		/* pending item exists?     */
	ipt = sspt->ss_stk[0];			/*   then use it	    */
    else {					/* nothing pending?	    */
	ipt = readitem(sspt->ss_str, NULL);	/*   read next item in      */
	sspt->ss_stk[0] = ipt;			/*   and save for later     */
    }
    return (ipt);				/* supply item to caller    */
}

local itemptr finditem(strstkptr sspt, string tag)
{
    itemptr sptr, *ivec;

    sptr = sspt->ss_stk[sspt->ss_stp];		/* get set from stack	    */
    ivec = (itemptr *) ItemDat(sptr);		/* get vect of items	    */
    while (*ivec != NULL) {			/* loop over set items      */
	if (streq(tag, ItemTag(*ivec)))		/*   found named item?      */
	    break;				/*     done with loop	    */
	ivec++;					/*   on to next item	    */
    }
    return (*ivec);				/* return item or NULL      */
}

/*
 * READITEM: read a simple or compound item.
 */

local itemptr readitem(stream str, itemptr first)
{
    itemptr ip, ibuf[MaxSetLen], *bufp, np, res;

    ip = first != NULL ? first : getitem(str);	/* use 1st or next item     */
    if (ip == NULL ||				/* EOF detected by getitem  */
	  ! streq(ItemTyp(ip), SetType))	/* or item not a set?       */
	return (ip);				/*   just return it	    */
    bufp = &ibuf[0];				/* prepare item buffer	    */
    for ( ; ; ) {				/* loop reading items in    */
	if (bufp >= &ibuf[MaxSetLen])		/*   no room for next?	    */
	    error("readitem: set %s: buffer overflow", ItemTag(ip));
	np = getitem(str);		        /*   look at next item	    */
	if (np == NULL)				/*   at end of file?	    */
	    error("readitem: set %s: unexpected EOF", ItemTag(ip));
	if (streq(ItemTyp(np), TesType))	/*   end of set item?	    */
	    break;				/*     quit input loop	    */
	*bufp++ = readitem(str, np);	  	/*   read next component    */
    }
    *bufp = NULL;				/* terminate item vector    */
    res = makeitem(scopy(SetType), scopy(ItemTag(ip)),
		   copxstr(ibuf, sizeof(itemptr)), NULL);
						/* construct compound item  */
    freeitem(ip, TRUE);				/* reclaim orig. header     */
    freeitem(np, TRUE);
    return (res);				/* return compound item     */
}

/*
 * GETITEM: read item from stream; return ptr to item or NULL.
 */

local itemptr getitem(stream str)
{
    itemptr ipt;

    ipt = gethdr(str);				/* try reading header in    */
    if (ipt == NULL)				/* did gethdr detect EOF?   */
        return (NULL);				/*   return NULL on EOF     */
    if (! streq(ItemTyp(ipt), SetType) &&	/* if item does not start   */
	  ! streq(ItemTyp(ipt), TesType))	/* or terminate a set       */
	getdat(ipt, str);			/*   try reading data in    */
    return (ipt);                               /* return resulting item */
}

/*
 * GETHDR: read a item header from a stream.
 */

local itemptr gethdr(stream str)
{
    short num;
    string typ, tag;
    int *dim, *ip;  /* ISSWAP */
    permanent bool firsttime = TRUE;

    if (fread(&num, sizeof(short), 1, str) != 1)/* read magic number*/
	return NULL;				/*   return NULL on EOF     */
    if (num == SingMagic || num == PlurMagic) {	/* new-style magic number?  */
	typ = (string) getxstr(str, sizeof(char));
						/*   read type string       */
	if (typ == NULL)			/*   check for EOF          */
	    error("gethdr: EOF reading type");
        swap = FALSE;
    }
#if defined(CHKSWAP)
    else {        /* ISSWAP */
        bswap((char *)&num,sizeof(short int),1);        /* swap the bytes */
        if (num == SingMagic || num == PlurMagic) {    /* test the swapped */
            if (firsttime)
                fprintf(stderr,"[filestruct: reading swapped]");
	    typ = (string) getxstr(str, sizeof(char));
						/*   read type string       */
	    if (typ == NULL)			/*   check for EOF          */
	        error("gethdr: EOF reading type");
            swap = TRUE;
            firsttime = FALSE;
        } else {
            bswap(&num,sizeof(short int),1);
            error("gethdr: bad magic: %o", num);
        }
    }
#else
    else {  
        error("gethdr: bad magic: %o",num);
    }
#endif
    if (! streq(typ, TesType)) {		/* is item tag next?        */
	tag = (string) getxstr(str, sizeof(char));
	if (tag == NULL)			/*   check for EOF          */
	    error("gethdr: EOF reading tag");
    } else
	tag = NULL;				/*   item is not tagged     */
    if (num == PlurMagic) {			/* are dimensions next?     */
	dim = (int *) getxstr(str, sizeof(int));
	if (dim == NULL)			/*   check for EOF          */
	    error("gethdr: EOF reading dimensions");
#if defined(CHKSWAP)
        if (swap) {     /* ISSWAP */
            ip = dim;
            while (*ip) {
                bswap((char *)ip,sizeof(int),1);
                ip++;
            }
        }
#endif
    } else
	dim = NULL;
    return (makeitem(typ, tag, NULL, dim));	/* return item less data    */
} /* gethdr */
/*
 * GETHDR: read a item header from a stream.
 */

bool qsf(stream str)
{
    short num;

    if (isatty(fileno(str)))   /* pipes and such cannot be structured files */
        return FALSE;

    if (fread(&num, sizeof(short), 1, str) != 1)/* read magic number        */
	return FALSE;				/*   return NULL on EOF     */
    if (num == SingMagic || num == PlurMagic) {	/* new-style magic number?  */
        return TRUE;
    }
#if defined(CHKSWAP)
    else {
        bswap(&num,sizeof(short int),1);        /* swap the bytes */
        if (num == SingMagic || num == PlurMagic) {    /* test the swapped */
            return TRUE;
        } else {
            return FALSE;
        }
    }
#else
    else {  
        return FALSE;
    }
#endif
} /* qsf */

/*
 * GETDAT: read item data from a stream.
 */

#define MaxReadNow  256

local void getdat(itemptr ipt, stream str)
{
    size_t dlen, elen;

    elen = eltcnt(ipt, 0);
    dlen = elen * ItemLen(ipt);                 /* count bytes of data	    */
#if 0
    if (dlen <= MaxReadNow) {			/* small enough to read?    */
#else
    if (dlen <= MaxReadNow || !strseek(str)) {  /* force read               */
#endif
	if(dlen<0) error("get_dat: dlen=%d",dlen);   
	ItemDat(ipt) = (byte *) calloc(dlen,1);	/*   then alloc space now   */
	if (ItemDat(ipt) == NULL)		/*   did alloc fail?	    */
	    error("getdat: no memory (%d bytes)", dlen);
	saferead(ItemDat(ipt), ItemLen(ipt), elen, str);
						/*   read data in now       */
    } else {					/* too big, so skip now     */
	ItemDat(ipt) = NULL;			/*   no data in core	    */
	ItemPos(ipt) = ftello(str);		/*   remember this place    */
	safeseek(str, dlen, 1);			/*   skip over data	    */
    }
} /* getdat */

/*
 * COPYFUN: select copy routine for given data types.
 */

local copyproc copyfun(string srctyp, string destyp)
{
    if (streq(srctyp, destyp))
	return copydata;
    if (streq(srctyp, FloatType) && streq(destyp, DoubleType))
	return (copyproc) copydata_f2d;
    if (streq(srctyp, DoubleType) && streq(destyp, FloatType))
	return (copyproc) copydata_d2f;
    return NULL;
} /* copyfun */

/*
 * COPYDATA - copy real or virtual data to assigned memory space
 */

local void copydata(
    void *vdat,
    int off,
    int len,
    itemptr ipt,
    stream str)
{
    char *src, *dat = (char *) vdat;
    off_t oldpos;
      
    off *= ItemLen(ipt);                        /* offset bytes from start  */
    if (ItemDat(ipt) != NULL) {			/* data already in core?    */
	src = (char *) ItemDat(ipt) + off;	/*   get pointer to source  */
    	len *= ItemLen(ipt);			/* number of bytes to copy  */
	while (--len >= 0)			/*   loop copying data      */
	    *dat++ = *src++;			/*     byte by byte         */
    } else {					/* time to read data in     */
	oldpos = ftello(str);                   /*   save current place     */
	safeseek(str, ItemPos(ipt) + off, 0);   /*   seek back to data      */
	saferead(dat, ItemLen(ipt), len, str);	/*   read data and check it */
	safeseek(str, oldpos, 0);               /*   reset file pointer     */
    }
} /* copydata */

local void copydata_f2d(
    double *dat,
    int off,
    int len,
    itemptr ipt,
    stream str)
{
    float *src;
    off_t oldpos;
      
    off *= ItemLen(ipt);
    if (ItemDat(ipt) != NULL) {			/* data already in core?    */
	src = (float *) ItemDat(ipt) + off;	/*   get pointer to source  */
	/*    	len *= ItemLen(ipt);		==> BUG */
	while (--len >= 0)			/*   loop converting data   */
	    *dat++ = (double) *src++;		/*     float to double      */
    } else {					/* time to read data in     */
	oldpos = ftello(str);                   /*   save this position     */
	safeseek(str, ItemPos(ipt) + off, 0);	/*   seek back to data      */
	while (--len >= 0)			/*   loop reading data      */
	    *dat++ = (double) getflt(str);	/*     float to double      */
	safeseek(str, oldpos, 0);               /*   reset file pointer     */
    }
} /* copydata_f2d */

local void copydata_d2f(
    float *dat,
    int off,
    int len,
    itemptr ipt,
    stream str)
{
    double *src;
    off_t oldpos;
      
    off *= ItemLen(ipt);
    if (ItemDat(ipt) != NULL) {			/* data already in core?    */
	src = (double *) ItemDat(ipt) + off;	/*   get pointer to source  */
    	/* len *= ItemLen(ipt);		BUG <===	*/
	while (--len >= 0)			/*   loop converting data   */
	    *dat++ = (float) *src++;		/*     double to float      */
    } else {					/* time to read data in     */
	oldpos = ftello(str);                   /*   save this position     */
	safeseek(str, ItemPos(ipt) + off, 0);	/*   seek back to data      */
	while (--len >= 0)			/*   loop reading data      */
	    *dat++ = (float) getdbl(str);	/*     double to float      */
	safeseek(str, oldpos, 0);               /*   reset file pointer     */
    }
} /* copydata_d2f */

local float getflt(stream str)
{
    float x;

    saferead((char *)&x, sizeof(float), 1, str);
    return x; 
}

local double getdbl(stream str)
{
    double x;

    saferead((char *)&x, sizeof(double), 1, str);
    return x;
}

local void saferead(
    void *dat,
    int siz,
    int cnt,
    stream str)
{
    if (fread(dat, siz, cnt, str) != cnt)
	error("saferead: error calling fread %d*%d bytes", siz, cnt);
#if defined(CHKSWAP)
    if (swap) bswap(dat,siz,cnt);
#endif
}

local void safeseek(
    stream str,
    off_t offset,
    int key)	      /* 0: start of file, 1: current point, 2: end */
{
    if (fseeko(str, offset, key) == -1)
	error("safeseek: error calling fseeko %d bytes from %d",
	      offset, key);
}

/************************************************************************/
/*                               UTILITIES                              */
/************************************************************************/

/*
 * ELTCNT: compute number of basic elements in subspace of item.
 */

local int eltcnt(
    itemptr ipt,            	/* pointer to item w/ possible vector dims */
    int skp)                	/* num of dims to skip, starting with dimN */
{
    register int prod, *ip;

    prod = 1;                                   /* scalers have one         */
    if (ItemDim(ipt) != NULL) {                 /* a vectorized item?       */
        for (ip = ItemDim(ipt); *ip != 0; ip++) /*   loop over dimensions   */
            if (--skp < 0)                      /*     past 1st skp dims?   */
                prod *= *ip;                    /*       include this dim   */
    }
    return (prod);				/* return product of dims   */
}

/*
 * DATLEN: compute length in bytes of subspace of item.
 */

local size_t datlen(
    itemptr ipt,            	/* pointer to item w/ possible vector dims */
    int skp)                	/* num of dims to skip, starting with dimN */
{
    return (ItemLen(ipt) * eltcnt(ipt, skp));	/* find byte-len needed     */
}

/*
 * MAKEITEM: construct an item from components.
 */

local itemptr makeitem(
    string typ,			/* data type of generated item */
    string tag,     	        /* tag for generated item */
    void *dat,              	/* pointer to data */
    int *dim)                 	/* dims, terminated with 0 */
{
    itemptr ipt;

    ipt = (itemptr) calloc(sizeof(item),1);	/* get space for item       */
    if (ipt == NULL)				/* check malloc worked      */
	error("makeitem: tag %s: malloc failed", tag);
    ItemTyp(ipt) = typ;				/* set type code string     */
    ItemLen(ipt) = baselen(typ);		/* set basic datum length   */
    ItemTag(ipt) = tag;				/* set item tag string      */
    if (dim != NULL && dim[0] != 0)             /* is item dimensioned?     */
        ItemDim(ipt) = dim;			/*   set dimension array    */
    else
	ItemDim(ipt) = NULL;			/*   clear out dimensions   */
    ItemDat(ipt) = dat;				/* set pointer to data      */
    ItemPos(ipt) = 0;				/* clear out file position  */
    return (ipt);                               /* return complete item     */
}

/*
 * FREEITEM: deallocate item and data.
 */
local void freeitem(
    itemptr ipt,			/* address of item to free */
    bool flg)			/* if true, free fields of item */
{
    itemptr *ivp;

    if (flg && ItemTyp(ipt) != NULL) {		/* does item have a type?   */
	if (streq(ItemTyp(ipt), SetType)) {	/*   free set recursively?  */
	    ivp = (itemptr *) ItemDat(ipt);	/*     get vector of items  */
	    if (ivp != NULL) 			/*     is data given?       */
		while (*ivp != NULL)		/*       loop over item set */
		    freeitem(*ivp++, TRUE);	/*         and free them up */
	}
	free(ItemTyp(ipt));			/*   free type string	    */
    }
    if (flg && ItemTag(ipt) != NULL)            /* is tag field set?        */
        free(ItemTag(ipt));                     /*   then free copy of tag  */
    if (flg && ItemDim(ipt) != NULL)
        free(ItemDim(ipt));
    if (flg && ItemDat(ipt) != NULL)
        free(ItemDat(ipt));
    free(ipt);                                  /* free item itself         */
}

/*
 * BASELEN: compute length of basic type in bytes.
 */

local typlen tl_tab[] = {
    { AnyType,	  sizeof(byte),   },
    { CharType,	  sizeof(char),   },
    { ByteType,	  sizeof(byte),   },
    { ShortType,  sizeof(short),  },
    { IntType,	  sizeof(int),    },
    { LongType,	  sizeof(long),   },
    { HalfpType,  sizeof(short),  },
    { FloatType,  sizeof(float),  },
    { DoubleType, sizeof(double), },
    { SetType,    0,              },
    { TesType,	  0,              },
    { NULL,	  0,              },
};

local int baselen(string typ)
{
    typlenptr tp;

    for (tp = tl_tab; tp->tl_typ != NULL; tp++)	/* loop over basic types    */
	if (streq(typ, tp->tl_typ))		/*   found type we want?    */
	    return (tp->tl_len);		/*     return byte length   */
    error("baselen: type %s unknown", typ);	/* bad type string          */
    return (0);                                 /* will never get here      */
}

/************************************************************************/
/*                             FINDING                                  */
/************************************************************************/
local string findtype(string *a, string type)
{
    int i;
    char *cp;

    for (i=0; a[i]; i++) {
        cp = a[i];
        if (streq(type,DoubleType) && *cp=='d')
            return cp;
        else if (streq(type,FloatType) && *cp=='f')
            return cp;
        else if (streq(type,HalfpType) && *cp=='h')
            return cp;
        else if (streq(type,IntType) && *cp=='i')
            return cp;
        else if (streq(type,IntType) && *cp=='s')
            return cp;
    }
    return NULL;
}

/************************************************************************/
/*                             STREAM TABLE                             */
/************************************************************************/

/*
 * FINDSTREAM: scan strtable for entry associated with stream.
 * If none is found, a new entry is initialized.
 * The table is small, no hash table is needed, although a quick
 * lookup is provided by first comparing it with the one used in 
 * the previous findstream() call.
 * A small degrading in performance may be the result in alternating I/O
 */

local strstk strtable[StrTabLen], *last = NULL;

local strstkptr findstream(stream str)
{
    strstkptr stfree, sspt;

    if (last!=NULL && last->ss_str == str)
        return(last);
    stfree = NULL;				/* remember free slot	    */
    for (sspt = strtable; sspt < strtable+StrTabLen; sspt++) {
						/* loop over the table	    */
	if (sspt->ss_str == str) {		/*   found that stream?	    */
            last = sspt;
	    return (sspt);			/*     then return slot	    */
        }
	if (stfree == NULL && sspt->ss_str == NULL)
						/*   first empty slot?	    */
	    stfree = sspt;			/*     save free slot ptr   */
    }
    if (stfree == NULL)				/* no free slot left?	    */
      error("findstream: no free slots, StrTabLen=%d",StrTabLen);

    stfree->ss_str = str;			/* init saved stream	    */
    stfree->ss_stk[0] = NULL;			/* clear pending item	    */
    stfree->ss_stp = -1;			/* empty item stack	    */
    stfree->ss_seek = TRUE;			/* permit seeks on stream   */
#if defined(RANDOM)
    stfree->ss_ran = NULL;                      /* mark as no item random   */
    stfree->ss_pos = 0L;                        /* set at start of file     */
#endif
    last = stfree;                              /* mark for quick access    */
    return (stfree);				/* return new slot	    */
}

local void ss_push(strstkptr sspt, itemptr ipt)
{
    if (sspt->ss_stp++ == SetStkLen)		/* check stack overflow	    */
	error("get_set: Too many nested items");
    sspt->ss_stk[sspt->ss_stp] = ipt;		/* put set on stack	    */
}

local void ss_pop(strstkptr sspt)
{
    if (sspt->ss_stp == -1)			/* check stack underflow    */
	error("ss_pop: stream stack underflow");
    sspt->ss_stp--;				/* bump stack pointer	    */
}

/************************************************************************/
/*			USER STREAM CONTROL FUNCTIONS			*/
/************************************************************************/

/*
 * STRCLOSE: remove stream from strtable, free associated items, and close.
 */

void strclose(stream str)
{
    strstkptr sspt;

    sspt = findstream(str);			/* lookup associated entry  */
    if (sspt->ss_stp != -1)			/* dont close if incomplete */
	error("strclose: not at top level");
    if (sspt->ss_stk[0] != NULL)		/* anything on the stack?   */
	freeitem(sspt->ss_stk[0], TRUE);	/*   free bottom item	    */
    sspt->ss_str = NULL;			/* remove from strtable	    */
    last = NULL;                                /* also removed quick access*/
    strdelete(str,FALSE);                       /* delete file if scratch   */
    fclose(str);				/* and close it up for sure */
}

