/*
 * FILESTRUCT.H: parameters and definitions for structured binary files.
 */

/*
 * Type codes.
 */

#define AnyType    000          /* anything at all */
#define CharType   001          /* printable chars */
#define ByteType   002          /* unprintable chars */
#define ShortType  003          /* short integers */
#define IntType    004          /* standard integers */
#define LongType   005          /* long integers */
#define FloatType  006          /* short floating */
#define DoubleType 007          /* long floating */
#define SetType    010          /* begin set of items */
#define TesType    011          /* end set of items */

/*
 * Function declerations.
 */

void   strclose(/* str */);

bool   get_tag_ok(/* str, tag */);

void   get_data(/* str, tag, typ, dat, dimN, ..., dim1, 0 */);
string get_string(/* str, tag */);
void   get_set(/* str, tag */);
void   get_tes(/* str */);

void   put_data(/* str, tag, typ, dat, dimN, ..., dim1, 0 */);
void   put_string(/* str, tag, msg */);
void   put_set(/* str, tag */);
void   put_tes(/* str */);

/*
 * For the benefit of applications programs which include this file,
 * define some standard names for various synonyms.
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

/*
 * MaxTagLen, MaxVecDim, MaxSetLen: arbitrary storage limits. These may be
 * increased as necessary without rendering old data files obsolete.
 */

#define MaxTagLen  65		/* max tag length, limited for simplicity */
#define MaxVecDim   9		/* max num of vec dim, limited for safety */
#define MaxSetLen  65		/* max num of components in compound item */
