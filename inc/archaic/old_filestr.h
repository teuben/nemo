/*
 * FILESTR.H: parameters and definitions for structured binary files.
 * A structured binary file may be viewed on three levels: as a stream
 * of (1) 8-bit bytes, (2) simple items, or (3) compound items.
 */

/*
 * IDES: encode flags, basic type and byte-length of basic types.
 * A short is given over to store the length so longish structures,
 * etc, may someday be included as basic types.
 */

typedef struct {
    byte  idesflag;		/* flags for input/output operations */
    byte  idestype;             /* basic data type: char, int, ... */
    short idestlen;             /* byte-length of basic dataum */
} ides;

#define TagFlag    001          /* item has a tag field */
#define VecFlag    002          /* item is a vector type */

#define AnyType    000          /* anything at all */
#define CharType   001          /* printable chars */
#define ByteType   002          /* unprintable chars */
#define ShortType  003          /* short integers */
#define IntType    004          /* standard integers */
#define LongType   005          /* long integers */
#define FloatType  006          /* short floating */
#define DoubleType 007          /* long floating */
#define BeginType  010          /* begin compound item */
#define EndType    011          /* end of compound item */

#define NBASETYPE (1 + EndType - AnyType)      /* number of basic types */

/*
 * ITEM: contents of an item, as stored internally.
 */

typedef struct {
    ides  itemdesc;             /* 4-byte flags, type, length */
    char *itemtag;              /* item ident, NULL-terminated */
    int  *itemdim;              /* vector dim's, 0-terminated */
    byte *itemdat;              /* address of actual data */
} item, *itemptr;

#define ItemDesc(ipt) ((ipt)->itemdesc)
#define ItemFlag(ipt) ((ipt)->itemdesc.idesflag)
#define ItemType(ipt) ((ipt)->itemdesc.idestype)
#define ItemTLen(ipt) ((ipt)->itemdesc.idestlen)
#define ItemTag(ipt)  ((ipt)->itemtag)
#define ItemDim(ipt)  ((ipt)->itemdim)
#define ItemDat(ipt)  ((ipt)->itemdat)

/*
 * MTAGLEN, MVECDIM, MSTRLEN: arbitrary storage limits. These may be
 * increased as necessary without rendering old data files obsolete.
 */

#define MTAGLEN  65             /* max tag length, limited for simplicity */
#define MVECDIM   9             /* max num of vec dim, limited for safety */
#define MSTRLEN	 65		/* max num of components in compound item */

/*
 * Function declarations.
 */

	/* high-level primatives */

bool writeitem(/* stream, itemptr */);

itemptr readitem(/* stream */);

itemptr finditem(/* stream, string */);

itemptr consitem(/* string, itemptr, . . ., NULL */);

itemptr scanitem(/* itemptr, string, . . ., NULL */);

void freeitem(/* itemptr, bool */);

	/* data-oriented extensions */

itemptr makeitem(/* byte, string, byte*, int, . . ., 0 */);

void putval(/* stream, byte, string, byte*, int, . . ., 0 */);

void getval(/* stream, byte, string, byte*, int, . . ., 0 */);

void putstr(/* stream, string, string */);

string getstr(/* stream, string */);

	/* low-level primatives */

bool putitem(/* stream, itemptr */);

itemptr getitem(/* stream */);

itemptr nxtitem(/* stream */);

	/* miscellaneous functions */

int eltcnt(/* itemptr, int */);

int datlen(/* itemptr, int */);

/*
 * For the benefit of applications programs which include
 * this file, define RealType depending on the precision,
 * and HeadlineTag to be the standard headline.
 */

#define BoolType ShortType

#ifndef SINGLEPREC
#  define RealType DoubleType
#else
#  define RealType FloatType
#endif

#ifndef HeadlineTag
#  define HeadlineTag "Headline"
#endif
