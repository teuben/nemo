/*
 * devices
 */

#ifndef VOGL
#define	VOGL
#endif

/*
 * valuator values
 */
#define MOUSEX		1
#define MOUSEY		2

/*
 * keyboard
 */
#define AKEY		'A'
#define BKEY		'B'
#define CKEY		'C'
#define DKEY		'D'
#define EKEY		'E'
#define FKEY		'F'
#define GKEY		'G'
#define HKEY		'H'
#define IKEY		'I'
#define JKEY		'J'
#define KKEY		'K'
#define LKEY		'L'
#define MKEY		'M'
#define NKEY		'N'
#define OKEY		'O'
#define PKEY		'P'
#define QKEY		'Q'
#define RKEY		'R'
#define SKEY		'S'
#define TKEY		'T'
#define UKEY		'U'
#define VKEY		'V'
#define WKEY		'W'
#define XKEY		'X'
#define YKEY		'Y'
#define ZKEY		'Z'
#define ZEROKEY		'0'
#define ONEKEY		'1'
#define TWOKEY		'2'
#define THREEKEY	'3'
#define FOURKEY		'4'
#define FIVEKEY		'5'
#define SIXKEY		'6'
#define SEVENKEY	'7'
#define EIGHTKEY	'8'
#define NINEKEY		'9'

#define SPACEKEY        ' '
#define SEMICOLONKEY    ';'
#define PERIODKEY       '.'
#define COMMAKEY        ','
#define QUOTEKEY        '\''
#define MINUSKEY        '-'
#define BACKSLASHKEY    '\\'
#define EQUALKEY        '='
#define LEFTBRACKETKEY  '['
#define RIGHTBRACKETKEY ']'

#define BACKSPACEKEY    '\010'
#define TABKEY          '\011'
#define LINEFEEDKEY     '\012'
#define RETKEY          '\015'
#define DELKEY          '\020'
#define ESCKEY		'\033'

#define	KEYBD		257

/*
 * mouse buttons
 */
#define MOUSE1		258
#define MOUSE2		259
#define MOUSE3		260
#define LEFTMOUSE	260
#define MIDDLEMOUSE	259
#define RIGHTMOUSE	258

/*
 * miscellany
 */
#define REDRAW		261
#define INPUTCHANGE	262

#define MAXDEV		262

/*
 * max size of device table.
 */
#define MAXDEVTABSIZE	34

extern long	qread();
extern Boolean	isqueued();

extern void	qdevice();
extern void	unqdevice();
extern void	qreset();

