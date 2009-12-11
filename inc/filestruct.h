/*
 * FILESTRUCT.H: user definitions for filestruct routines.
 *	4-mar-91  added copy_item_cvt		PJT
 *	6-may-92  random data access routines	PJT
 *     12-apr-95  prototypes without ARGS       PJT
 *      2-jun-05  blocked I/O as a flavor of random I/O     PJT
 *     11-dec-09  half precision type                       PJT
 */
#ifndef _filestruct_h
#define _filestruct_h
/*
 * New-style basic item types.
 */

#define AnyType    "a"       /* anything at all */
#define CharType   "c"       /* printable chars */
#define ByteType   "b"       /* unprintable chars */
#define ShortType  "s"       /* short integers */
#define IntType    "i"       /* standard integers */
#define LongType   "l"       /* long integers */
#define HalfpType  "h"       /* half precision floating */
#define FloatType  "f"       /* short floating */
#define DoubleType "d"       /* long floating */
#define SetType    "("       /* begin compound item */
#define TesType    ")"       /* end of compound item */
/* Experimental Kludge */
#define StoryType  "["       /* begin of a story item (see starlab) */
#define YrotsType  "]"       /* end of a story item (see starlab) */
/*
 * MaxTagLen, MaxVecDim, MaxSetLen: arbitrary storage limits. These may be
 * increased as necessary without rendering old data files obsolete.
 */

#define MaxTagLen  65		/* max tag length, limited for simplicity */
#define MaxVecDim   9		/* max num of vec dim, limited for safety */
#define MaxSetLen  65		/* max num of components in compound item */

/*
 * For the benefit of applications programs which include this file,
 * define standard synonyms for bools and reals (see stdinc.h).
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
 * Publically accessible routines.
 */

extern void copy_item (stream, stream, string);
extern void copy_item_cvt (stream, stream, string, string *);

extern void put_set (stream, string);
extern void put_tes (stream, string);
extern void put_string ( stream, string , string );
extern void put_data ( stream, string, string, void *, int, ...);
extern void put_data_sub ( stream, string, string, void *, int *, bool);

extern void get_set ( stream str, string tag );
extern void get_tes ( stream str, string tag );
extern string get_string ( stream str, string tag );
extern void get_data ( stream, string, string, void *, int, ...);

extern void get_data_coerced ( stream, string, string, void *, int, ...);
			 
extern void get_data_sub ( stream, string, string, void *, int *, bool);
		     
extern bool get_tag_ok ( stream, string);
extern bool skip_item ( stream);
extern string *list_tags ( stream);
extern string get_type ( stream, string);
extern int *get_dims ( stream, string);
extern size_t get_dlen ( stream, string);

extern void strclose ( stream);

extern void get_data_set     ( stream , string , string , int,  ...);
extern void get_data_tes     ( stream , string  );
extern void get_data_ran     ( stream , string , void *, int , int );
extern void get_data_blocked ( stream , string , void *, int);

extern void put_data_set     ( stream , string , string , int,  ...);
extern void put_data_tes     ( stream , string );
extern void put_data_ran     ( stream , string , void *, int , int  );
extern void put_data_blocked ( stream , string , void *, int );

extern bool qsf ( stream );
#endif
