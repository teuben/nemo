/*
 * FILESECRET.H: parameters and definitions for secret (user-invisible)
 *		 part of structured file package.
 * V 1.0: Joshua Barnes 4/87	basic I/O operators implemented,
 * V 1.X: Lyman Hurd 8/87	added f-d coercion, deferred input,
 * V 2.0: Joshua Barnes 4/88	new types, operators, external format.
 *	  4-oct-90		bug in ss_stk[SetStkLen] removed: 1 more
 * V 2.1: Peter Teuben 13-oct-90  random access
 * V 2.2: 26-feb-94   ansi				PJT
 *	  14-mar-95   minor cleanup
 */
 
#define RANDOM  /* allow random access */
#define CHKSWAP /* allow mixed endian datasets - this is a bit dangerous */

/*
 * New-style magic numbers, for (bigendian) FITS type machines (like SUN)
 * Note: If they appear swapped, filestruct() runs in swapped mode !
 */

#define SingMagic  ((011<<8) + 0222)		/* singular items */
#define PlurMagic  ((013<<8) + 0222)		/* plural items */

/*
 * ITEM: structure representing data-token.
 */

typedef struct {
    string itemtyp;		/* type string, listed in filestruct.h */
    long   itemlen;		/* length associated with above type */
    string itemtag;		/* name given this item by application */
    int   *itemdim;		/* int-string of dimensions, or NULL */
    void  *itemdat;		/* the real goodies, if any, or NULL */
    off_t  itempos;		/* where the item began in stream (i/o) */
} item, *itemptr;    

#define ItemTyp(ip)  ((ip)->itemtyp)
#define ItemLen(ip)  ((ip)->itemlen)
#define ItemTag(ip)  ((ip)->itemtag)
#define ItemDim(ip)  ((ip)->itemdim)
#define ItemDat(ip)  ((ip)->itemdat)
#define ItemPos(ip)  ((ip)->itempos)

/*
 * STRSTK: structure used to associate stream with item stack.
 */

#define SetStkLen  8
#define StrTabLen 16

typedef struct {
    stream  ss_str;		    /* pointer to stdio stream */
    itemptr ss_stk[SetStkLen+1];    /* stack of assoc. items (extra dummy) */
    int     ss_stp;		    /* item stack pointer */
    bool    ss_seek;		    /* permit seeks on this stream ? */
#if defined(RANDOM)
    off_t   ss_pos;                 /* tail of file, in case random access */
    itemptr ss_ran;                 /* pointer to random access item */
#endif
} strstk, *strstkptr;

/*
 * TYPLEN: structure used to associate length with item type.
 */

typedef struct {
    string tl_typ;			/* type string (see filestruct.h)   */
    int    tl_len;			/* length of elements in bytes      */
} typlen, *typlenptr;

/*
 * Function declarations.
 */

/* Local functions */

local bool writeitem   ( stream str, itemptr ipt );
local bool putitem     ( stream str, itemptr ipt );
local bool puthdr      ( stream str, itemptr ipt );
local bool putdat      ( stream str, itemptr ipt );
local itemptr scantag  ( strstkptr sspt, string tag );
local itemptr nextitem ( strstkptr sspt );
local itemptr finditem ( strstkptr sspt, string tag );
local itemptr readitem ( stream str, itemptr first );
local itemptr getitem  ( stream str );
local itemptr gethdr   ( stream str );
local void getdat      ( itemptr ipt, stream str );
local proc copyfun     ( string srctyp, string destyp );
local void copydata    ( void *dat, int, int, itemptr ipt, stream str );
local void copydata_f2d ( double *dat, int, int, itemptr ipt, stream str );
local void copydata_d2f ( float  *dat, int, int, itemptr ipt, stream str );
local float getflt     ( stream str );
local double getdbl    ( stream str );
local void saferead    ( void *dat, int siz, int cnt, stream str );
local void safeseek    ( stream str, off_t offset, int key );
local int eltcnt       ( itemptr ipt, int skp );
local size_t datlen    ( itemptr ipt, int skp );
local itemptr makeitem ( string typ, string tag, void *dat, int *dim );
local void freeitem    ( itemptr ipt, bool flg);
local int baselen      ( string typ );
local strstkptr findstream ( stream str );
local void ss_push     ( strstkptr sspt, itemptr ipt );
local void ss_pop      ( strstkptr sspt );
local string findtype  ( string *a, string type );


#if defined(CHKSWAP)
 local bool swap=FALSE;
#endif

