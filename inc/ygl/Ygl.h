#ifndef	_YGL_INCLUDED_
#define	_YGL_INCLUDED_

#include <sys/types.h>

#define	BLACK			0
#define	WHITE			1
#define	GREEN			2
#define	YELLOW			3
#define	BLUE			4
#define	MAGENTA			5
#define	CYAN			6
#define	RED			7

#define LEFTMOUSE		103
#define MIDDLEMOUSE		102
#define RIGHTMOUSE		101
#define MENUBUTTON		RIGHTMOUSE

#define MOUSEX			266
#define MOUSEY			267

#define KEYBD			513
#define REDRAW			528
#define	INPUTCHANGE		534
#define WINCLOSE		537

typedef char			Int8;
typedef	unsigned char		Uint8;
typedef short			Int16;
typedef	unsigned short		Uint16;
typedef long			Int32;
typedef	unsigned long		Uint32;
typedef	float			Float32;
typedef	double			Float64;
typedef	char			Char8;
typedef char			Void;
typedef	Int16			Angle;
typedef	Uint8			Byte;
typedef	Uint16			Colorindex;
typedef	Float32			Coord;
typedef	Uint16			Device;
typedef	Int32			Icoord;
typedef	Uint8			RGBvalue;
typedef	Int16			Scoord;
typedef	Int16			Screencoord;

extern void  arc 		(  Coord,  Coord,  Coord, Angle, Angle );
extern void  arci		( Icoord, Icoord, Icoord, Angle, Angle );
extern void  arcs		( Scoord, Scoord, Scoord, Angle, Angle );

extern void  arcf 		(  Coord,  Coord,  Coord, Angle, Angle );
extern void  arcfi		( Icoord, Icoord, Icoord, Angle, Angle );
extern void  arcfs		( Scoord, Scoord, Scoord, Angle, Angle );

extern void  charstr		( char * );

extern void  circ 		(  Coord,  Coord,  Coord );
extern void  circi		( Icoord, Icoord, Icoord );
extern void  circs		( Scoord, Scoord, Scoord );

extern void  circf 		(  Coord,  Coord,  Coord );
extern void  circfi		( Icoord, Icoord, Icoord );
extern void  circfs		( Scoord, Scoord, Scoord );

extern void  clear		( void );

extern void  cmov2 		(  Coord,  Coord );
extern void  cmov2i		( Icoord, Icoord );
extern void  cmov2s		( Scoord, Scoord );

extern void  cmode		( void );

extern void  color		( Colorindex );
extern void  concave		( long );
extern void  doublebuffer	( void );

extern void  draw2 		(  Coord,  Coord );
extern void  draw2i		( Icoord, Icoord );
extern void  draw2s		( Scoord, Scoord );

extern void  font		( Int16	);
extern void  gconfig		( void );
extern long  getbutton		( Device );
extern Int32 getcolor		( void );
extern Int32 getdescender	( void );
extern Int32 getheight		( void );
extern void  getmcolor		( Colorindex, Int16 *, Int16 *,	Int16 *	);
extern void  getorigin		( Int32	*, Int32 * );
extern Int32 getplanes		( void );
extern void  getsize		( Int32	*, Int32 * );
extern Int32 getvaluator	( Device );
extern void  gexit		( void );
extern void  ginit		( void );
extern void  gRGBcolor		( Int16 *, Int16 *, Int16 * );
extern Int32 gversion		( Char8[12] );
extern void  keepaspect		( Int32, Int32 );
extern void  loadXfont		( Int32	id, char * n );
extern void  mapcolor		( Colorindex, Int16, Int16, Int16 );
extern void  maxsize		( Int32, Int32 );
extern void  minsize		( Int32, Int32 );

extern void  move2 		(  Coord,  Coord );
extern void  move2i		( Icoord, Icoord );
extern void  move2s		( Scoord, Scoord );

extern void  pclos		( void );

extern void  pdr2 		(  Coord,  Coord );
extern void  pdr2i		( Icoord, Icoord );
extern void  pdr2s		( Scoord, Scoord );

extern void  pmv2 		(  Coord,  Coord );
extern void  pmv2i		( Icoord, Icoord );
extern void  pmv2s		( Scoord, Scoord );

extern void  pnt2 		(  Coord,  Coord );
extern void  pnt2i		( Icoord, Icoord );
extern void  pnt2s		( Scoord, Scoord );

extern void  polf2 		( Int32,  Coord	[][2] );
extern void  polf2i		( Int32, Icoord	[][2] );
extern void  polf2s		( Int32, Scoord	[][2] );

extern void  poly2 		( Int32,  Coord	[][2] );
extern void  poly2i		( Int32, Icoord	[][2] );
extern void  poly2s		( Int32, Scoord	[][2] );

extern void  prefposition	( Int32, Int32,	Int32, Int32 );
extern void  prefsize		( Int32, Int32 );
extern void  qdevice		( Device );
extern void  qenter		( Int16, Int16 );
extern Int32 qread		( Int16	* );
extern void  qreset		( void );
extern Int32 qtest		( void );

extern void  rect 		(  Coord,  Coord,  Coord,  Coord );
extern void  recti		( Icoord, Icoord, Icoord, Icoord );
extern void  rects		( Scoord, Scoord, Scoord, Scoord );

extern void  rectf 		(  Coord,  Coord,  Coord,  Coord );
extern void  rectfi		( Icoord, Icoord, Icoord, Icoord );
extern void  rectfs		( Scoord, Scoord, Scoord, Scoord );

extern void  reshapeviewport	( void );

extern void  rdr2 		(  Coord,  Coord );
extern void  rdr2i		( Icoord, Icoord );
extern void  rdr2s		( Scoord, Scoord );

extern void  RGBcolor		( Int16, Int16, Int16 );
extern void  RGBmode		( void );

extern void  rpdr2 		(  Coord,  Coord );
extern void  rpdr2i		( Icoord, Icoord );
extern void  rpdr2s		( Scoord, Scoord );

extern void  rmv2 		(  Coord,  Coord );
extern void  rmv2i		( Icoord, Icoord );
extern void  rmv2s		( Scoord, Scoord );

extern void  sbox 		(  Coord,  Coord,  Coord,  Coord );
extern void  sboxi		( Icoord, Icoord, Icoord, Icoord );
extern void  sboxs		( Scoord, Scoord, Scoord, Scoord );

extern void  sboxf 		(  Coord,  Coord,  Coord,  Coord );
extern void  sboxfi		( Icoord, Icoord, Icoord, Icoord );
extern void  sboxfs		( Scoord, Scoord, Scoord, Scoord );

extern void  singlebuffer	( void );
extern void  stepunit		( Int32, Int32 );
extern Int32 strwidth		( Char8	* );
extern void  swapbuffers	( void );
extern void  unqdevice		( Device );
extern void  winclose		( Int32	);
extern Int32 winget		( void );
extern Int32 winopen		( Char8 * );
extern void  winset		( Int32	);
extern void  wintitle		( Char8 * );
extern void  winmove		( Int32, Int32 );
extern void  winposition	( Int32, Int32, Int32, Int32 );

/* Extensions: Routines from X not in gl. Mostly by MiSt (michael@hal6000.thp.Uni-Duisburg.DE) */

extern void  arcx 		(  Coord,  Coord,  Coord,  Coord, Angle, Angle );
extern void  arcxi		( Icoord, Icoord, Icoord, Icoord, Angle, Angle );
extern void  arcxs		( Scoord, Scoord, Scoord, Scoord, Angle, Angle );

extern void  arcxf 		(  Coord,  Coord,  Coord,  Coord, Angle, Angle );
extern void  arcxfi		( Icoord, Icoord, Icoord, Icoord, Angle, Angle );
extern void  arcxfs		( Scoord, Scoord, Scoord, Scoord, Angle, Angle );

#endif /* _YGL_INCLUDED_ */
