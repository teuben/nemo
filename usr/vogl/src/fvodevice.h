
c
c valuator values
c
	integer MOUSEX, MOUSEY
	parameter(MOUSEX = 1, MOUSEY = 2)
c
c keyboard and mouse buttons
c 
	integer AKEY, BKEY, CKEY, DKEY, EKEY, FKEY, GKEY, HKEY
	integer IKEY, JKEY, KKEY, LKEY, MKEY, NKEY, OKEY, PKEY
	integer QKEY, RKEY, SKEY, TKEY, UKEY, VKEY, WKEY, XKEY
	integer YKEY, ZKEY

	integer ZEROKE, ONEKEY, TWOKEY, THREEK, FOURKE, FIVEKE
	integer SIXKEY, SEVENK, EIGHTK, NINEKE 

	integer SPACEK, SEMICO, PERIOD, COMMAN, QUOTEK, MINUSK
	integer BACKSL, EQUALK, LEFTBR, RIGHTB

	integer BACKSP, TABKEY, LINEFE, RETKEY, DELKEY, ESCKEY

	integer KEYBD, MOUSE1, MOUSE2, MOUSE3, LEFTMO, RIGHTM
	integer MIDDLE, MAXDEV

	parameter(AKEY = 65, BKEY = 66, CKEY = 67, DKEY = 68)
	parameter(EKEY = 69, FKEY = 70, GKEY = 71, HKEY = 72)
	parameter(IKEY = 73, JKEY = 74, KKEY = 75, LKEY = 76)
	parameter(MKEY = 77, NKEY = 78, OKEY = 79, PKEY = 80)
	parameter(QKEY = 81, RKEY = 82, SKEY = 83, TKEY = 84)
	parameter(UKEY = 85, VKEY = 86, WKEY = 87, XKEY = 88)
	parameter(YKEY = 89, ZKEY = 90)

	parameter(ZEROKE = 48, ONEKEY = 49, TWOKEY = 50, THREEK = 51)
	parameter(FOURKE = 52, FIVEKE = 53, SIXKEY = 54, SEVENK = 55)
	parameter(EIGHTK = 56, NINEKE = 57)

	parameter(SPACEK = 32, SEMICO = 59, PERIOD = 46, COMMAN = 44)
	parameter(QUOTEK = 39, MINUSK = 45, BACKSL = 92, EQUALK = 61)
	parameter(LEFTBR = 91, RIGHTB = 93)

	parameter(BACKSP = 8, TABKEY = 9, LINEFE = 10, RETKEY = 13)
	parameter(DELKEY = 127, ESCKEY = 27)

	parameter(KEYBD = 257, MOUSE1 = 258, MOUSE2 = 259, MOUSE3 = 260)
	parameter(LEFTMO = 260, MIDDLE = 259, RIGHTM = 258)
c
c miscellany
c
	integer REDRAW, INPUTC
	parameter(REDRAW = 261, INPUTC = 262)

	parameter(MAXDEV = 260)

	integer getbut, getval, getpla, qread
	logical qtest
	external getbut, getval, getpla, qread, qtest
