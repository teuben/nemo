
#ifdef PC	/* Stupid pox head crap */
char	*vallocate();
char	*malloc();
#endif

/*
 * VOGL is always defined if a header file is from the 
 * VOGL library. In cases where you do use some VOGL
 * initialisation routines like vinit, just put #ifdef VOGL...
 * around.
 */
#ifndef VOGL
#define	VOGL
#endif

#ifndef TRUE
#define	TRUE	1
#endif

#ifndef FALSE
#define	FALSE	0
#endif

/*
 * Misc defines...
 */
#define	FLAT	0
#define SMOOTH	1
#define GD_XPMAX 1
#define GD_YPMAX 2

/*
 * standard colour indices
 */
#define	BLACK		0
#define	RED		1
#define	GREEN		2
#define	YELLOW		3
#define	BLUE		4
#define	MAGENTA		5
#define	CYAN		6
#define	WHITE		7

/*
 * when (if ever) we need the precision
 */
#ifdef DOUBLE
#define	float	double
#endif

/*
 * How to convert degrees to radians
 */
#define	PI	3.14159265358979
#define D2R	(PI / 180.0)

/*
 * miscellaneous typedefs and type defines
 */
typedef float	Vector[4];
typedef float	Matrix[4][4];
typedef float	Tensor[4][4][4];
typedef short	Angle;
typedef float	Coord;
typedef long	Icoord;
typedef short	Scoord;
typedef long	Object;
typedef short	Screencoord;
typedef long	Boolean;
typedef unsigned short	Linestyle;

typedef unsigned short	Device;

typedef unsigned short	Colorindex;


/*
 * when register variables get us into trouble
 */
#ifdef NOREGISTER
#define	register
#endif

/*
 * max number of vertices in a ploygon
 */
#define	MAXVERTS	128

/*
 * object definitions
 */
#define MAXENTS		101		/* size of object table */
#define	MAXTOKS		100		/* num. of tokens alloced at once in
					   an object  */

/*
 * Polygon fill modes for "polymode"
 */
#define PYM_POINT	0
#define PYM_LINE	0
#define PYM_FILL	1
#define PYM_HOLLOW	1

/*
 * functions which can appear in objects
 */
#define	ARC		1
#define	CALLOBJ		3
#define	CIRCLE		5
#define	CLEAR		6
#define	COLOR		7
#define	DRAW		8
#define	DRAWSTR		10
#define	VFONT		12
#define	LOADMATRIX	15
#define	MAPCOLOR	16
#define	MOVE		17
#define	MULTMATRIX	18
#define	POLY		19
#define	POPATTRIBUTES	22
#define	POPMATRIX	23
#define	POPVIEWPORT	24
#define	PUSHATTRIBUTES	25
#define	PUSHMATRIX	26
#define	PUSHVIEWPORT	27
#define	RCURVE		28
#define	RPATCH		29
#define	SECTOR		30
#define	VIEWPORT	33
#define	BACKBUFFER	34
#define	FRONTBUFFER	35
#define	SWAPBUFFERS	36
#define	BACKFACING	37
#define	TRANSLATE	38
#define	ROTATE		39
#define	SCALE		40

#define	ARCF		41
#define	CIRCF		42
#define	POLYF		43
#define	RECTF		44
#define	POLYMODE	45
#define	CMOV		46
#define	LINESTYLE	47
#define	LINEWIDTH	48

/*
 * Non standard call...
 */
#define	VFLUSH		70

/*
 * States for bgn* and end* calls
 */
#define	NONE		0	/* Just set current spot */
#define	VPNT		1	/* Draw dots		 */
#define	VLINE		2	/* Draw lines		 */
#define	VCLINE		3	/* Draw closed lines	 */
#define	VPOLY		4	/* Draw a polygon 	 */
#define VTMESH		5       /* Draw a triangular mesh*/
#define VQSTRIP		6       /* Draw a quadralateral mesh*/

/*
 * data types for object tokens
 */
typedef union tk {
	int		i;
	float		f;
} Token;

typedef struct tls {
	int		count;
	Token		*toks;
	struct tls	*next;
} TokList;

/*
 * double buffering modes.
 */
#define	SINGLE		1

/*
 * attributes
 */
typedef struct {
	char		backface,
			mode;			/* Which mode are we in */
	int		color;
	int		fontnum;
	Linestyle	ls;			
	short		lw;			/* Linewidth */
} Attribute;

/*
 * viewport
 */
typedef struct vp {
	float	left;
	float	right;
	float	bottom;
	float	top;
} Viewport; 

/*
 * stacks
 */
typedef	struct	ms {	/* Matrix stack entries	*/
	Matrix		m;
	struct	ms	*back;
} Mstack;

typedef	struct	as {	/* Attribute stack entries */
	Attribute	a;
	struct	as	*back;
} Astack;

typedef	struct	vs {	/* Viewport stack entries */
	Viewport	v;
	struct	vs	*back;
} Vstack;

/*
 * vogle device structures
 */
typedef struct dev {
	char	*devname;		/* name of device */
	char	*large,			/* name of large font */
		*small;			/* name of small font */
	int	(*Vbackb)(),		/* Set drawing in back buffer */
		(*Vchar)(),		/* Draw a hardware character */
		(*Vcheckkey)(),		/* Ckeck if a key was hit */
		(*Vclear)(),		/* Clear the screen to current color */
		(*Vcolor)(),		/* Set current color */
		(*Vdraw)(),		/* Draw a line */
		(*Vexit)(),		/* Exit graphics */
		(*Vfill)(),		/* Fill a polygon */
		(*Vfont)(),		/* Set hardware font */
		(*Vfrontb)(),		/* Set drawing in front buffer */
		(*Vgetkey)(),		/* Wait for and get the next key hit */
		(*Vinit)(),		/* Initialise the device */
		(*Vlocator)(),		/* Get mouse/cross hair position */
		(*Vmapcolor)(),		/* Set color indicies */
		(*Vsetls)(),		/* Set linestyle */
		(*Vsetlw)(),		/* Set linewidth */
		(*Vstring)(),		/* Draw a hardware string */
		(*Vswapb)(),		/* Swap front and back buffers */
		(*Vsync)();		/* Sync display */
} DevEntry;

typedef struct vdev {
	char		initialised,
			clipoff,
			inobject,
			inpolygon,
			fill,			/* polygon filling */
			cpVvalid,		/* is the current device position valid */
			sync,			/* Do we syncronise the display */
			inbackbuffer,		/* are we in the backbuffer */
			clipplanes;		/* active clipping planes */
	void		(*pmove)(),		/* Polygon moves */
			(*pdraw)();		/* Polygon draws */
	TokList		*tokens;		/* ptr to list of tokens for current object */
	Mstack		*transmat;		/* top of transformation stack */
	Astack		*attr;			/* top of attribute stack */
	Vstack		*viewport;		/* top of viewport stack */
	float		hheight, hwidth;	/* hardware character height, width */
	Vector		cpW,			/* current postion in world coords */
			cpWtrans,		/* current world coords transformed */
			upvector;		/* world up */
	int		depth,			/* # bit planes on screen */
			maxVx, minVx,
			maxVy, minVy,
			sizeX, sizeY, 		/* size of square on screen */
			sizeSx, sizeSy,		/* side in x, side in y (# pixels) */
			cpVx, cpVy;
	DevEntry	dev;
	float		savex,			/* Where we started for v*() */
			savey,
			savez;
	char		bgnmode;		/* What to do with v*() calls */
	int		save;			/* Do we save 1st v*() point */

	char		*wintitle;		/* window title */

	char		*devname;		/* pointer to device name */

	Matrix		tbasis, ubasis, *bases; /* Patch stuff */
	
	char		*enabled;		/* pointer to enabled devices mask */
	int		maxfontnum;

	char		alreadyread;		/* queue device stuff */
	char		kbdmode;		/* are we in keyboard mode */
	char		mouseevents;		/* are mouse events enabled */
	char		kbdevents;		/* are kbd events enabled */
	int		devno, data;

	int		concave;		/* concave polygons? */
} VDevice;

extern VDevice	vdevice;		/* device structure */

#define	V_X	0			/* x axis in cpW */
#define	V_Y	1			/* y axis in cpW */
#define	V_Z	2			/* z axis in cpW */
#define	V_W	3			/* w axis in cpW */

/*
 * function definitions
 */

/*
 * arc routines
 */
extern void	arcprecision();
extern void	circleprecision();
extern void	arc();
extern void	arcs();
extern void	arci();
extern void	arcf();
extern void	arcfs();
extern void	arcfi();
extern void	circ();
extern void	circs();
extern void	circi();
extern void	circf();
extern void	circfs();
extern void	circfi();

/*
 * attr routines
 */
extern void	popattributes();
extern void	pushattributes();

/*
 * curve routines
 */
extern void	curvebasis();
extern void	curveprecision();
extern void	rcrv();
extern void	crv();
extern void	crvn();
extern void	rcrvn();
extern void	curveit();

/*
 * draw routines
 */
extern void	draw();
extern void	draws();
extern void	drawi();
extern void	draw2();
extern void	draw2s();
extern void	draw2i();
extern void	rdr();
extern void	rdrs();
extern void	rdri();
extern void	rdr2();
extern void	rdr2s();
extern void	rdr2i();
extern void	bgnline();
extern void	endline();
extern void	bgnclosedline();
extern void	endclosedline();

/*
 * device routines
 */
extern void	qdevice();
extern void	unqdevice();
extern long	qread();
extern void	qreset();
extern long	qtest();
extern Boolean	isqueued();

extern void	gexit();
extern void	gconfig();
extern void	shademodel();
extern long	getgdesc();
extern long	winopen();
extern void	ginit();
extern void	gconfig();
extern long	getvaluator();
extern Boolean	getbutton();
extern void	getdev();
extern void	clear();
extern void	colorf();
extern void	color();
extern void	mapcolor();
extern long	getplanes();

extern void	vinit();
extern void	voutput();
extern void	verror();
extern void	vnewdev();
extern char	*vgetdev();

/*
 * mapping routines
 */
extern int	WtoVx();
extern int	WtoVy();
extern void	CalcW2Vcoeffs();

/*
 * general matrix and vector routines
 */
extern void	mult4x4();
extern void	copymatrix();
extern void	identmatrix();
extern void	copytranspose();

extern void	multvector();
extern void	copyvector();
extern void	premultvector();

/*
 * matrix stack routines
 */
extern void	getmatrix();
extern void	popmatrix();
extern void	loadmatrix();
extern void	pushmatrix();
extern void	multmatrix();

/*
 * move routines
 */
extern void	move();
extern void	moves();
extern void	movei();
extern void	move2();
extern void	move2s();
extern void	move2i();
extern void	rmv();
extern void	rmvs();
extern void	rmvi();
extern void	rmv2();
extern void	rmv2s();
extern void	rmv2i();

/*
 * object routines
 */
extern Boolean	isobj();
extern long	genobj();
extern void	delobj();
extern void	makeobj();
extern void	callobj();
extern void	closeobj();
extern long	getopenobj();
extern Token	*newtokens();

/*
 * patch routines.
 */
extern void	defbasis();
extern void	patchbasis();
extern void	patchcurves();
extern void	patchprecision();
extern void	patch();
extern void	rpatch();

/*
 * point routines
 */
extern void	pnt();
extern void	pnts();
extern void	pnti();
extern void	pnt2();
extern void	pnt2s();
extern void	pnt2i();
extern void	bgnpoint();
extern void	endpoint();

/*
 * v routines
 */
extern void	v4f();
extern void	v3f();
extern void	v2f();
extern void	v4d();
extern void	v3d();
extern void	v2d();
extern void	v4i();
extern void	v3i();
extern void	v2i();
extern void	v4s();
extern void	v3s();
extern void	v2s();

/*
 * polygon routines.
 */
extern void	concave();
extern void	backface();
extern void	frontface();
extern void	polymode();
extern void	poly2();
extern void	poly2i();
extern void	poly2s();
extern void	polyi();
extern void	polys();
extern void	polf2();
extern void	polf2i();
extern void	polf2s();
extern void	polfi();
extern void	polfs();
extern void	poly();
extern void	polf();
extern void	pmv();
extern void	pmvi();
extern void	pmv2i();
extern void	pmvs();
extern void	pmv2s();
extern void	pmv2();
extern void	pdr();
extern void	rpdr();
extern void	rpdr2();
extern void	rpdri();
extern void	rpdr2i();
extern void	rpdrs();
extern void	rpdr2s();
extern void	rpmv();
extern void	rpmv2();
extern void	rpmvi();
extern void	rpmv2i();
extern void	rpmvs();
extern void	rpmv2s();
extern void	pdri();
extern void	pdr2i();
extern void	pdrs();
extern void	pdr2s();
extern void	pdr2();
extern void	pclos();
extern void	bgnpolygon();
extern void	endpolygon();

/*
 * rectangle routines
 */
extern void	rect();
extern void	recti();
extern void	rects();
extern void	rectf();
extern void	rectfi();
extern void	rectfs();

/*
 * tensor routines
 */
extern void multtensor();
extern void copytensor();
extern void premulttensor();
extern void copytensortrans();

/*
 * text routines
 */
extern void	font();
extern void	charstr();
extern void	cmov();
extern void	cmov2();
extern void	cmovi();
extern void	cmovs();
extern void	cmov2i();
extern void	cmov2s();
#ifdef OLD_GL
extern long	getwidth();
#endif
extern long	getheight();
extern long	strwidth();
extern void	getcpos();

/*
 * transformation routines
 */
extern void	scale();
extern void	translate();
extern void	rotate();
extern void	rot();

/*
 * window definition routines
 */
extern void	ortho();
extern void	ortho2();
extern void	lookat();
extern void	window();
extern void	polarview();
extern void	perspective();

/*
 * routines for manipulating the viewport
 */
extern void	viewport();
extern void	popviewport();
extern void	pushviewport();

/*
 * routines for retrieving the graphics position
 */
extern void	getgp();
extern void	getgpos();

/*
 * routines for handling the buffering
 */
extern void	backbuffer();
extern void	frontbuffer();
extern void	swapbuffers();
extern void	doublebuffer();

/*
 * routines for window sizing and positioning
 */
extern void	prefsize();
extern void	prefposition();

/*
 * Misc control routines
 */
extern void	vsetflush();
extern void	vflush();
