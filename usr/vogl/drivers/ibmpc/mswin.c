/****************************************************************************

    PROGRAM: MSWIN.C

    PURPOSE: MS-Windows driver for VOGLE/VOGL

    FUNCTIONS:

	WinMain() - calls initialization function, calls user's VOGLE program (main)
	InitApplication() - initializes window data and registers window
	InitInstance() - saves instance handle
	MainWndProc() - processes messages
	About() - processes messages for "About" dialog box
	mswin_msg_loop() - grabs and dispatches windows messages (until no more pending)
	mswin_cleanup() - frees objects created along the way in preparation for exit
	mswin_...() - VOGLE driver routines called by VOGLE
	swap(() - swaps two integer values

    COMMENTS:

	The GENERIC sample program from Windows SDK was used as the basis
	for this driver.

****************************************************************************/

#include "windows.h"				/* required for all Windows applications */
#include <stdio.h>
#include "mswin.h"				/* specific to this program */
#ifdef VOGLE
#include "vogle.h"
#else
#include "vogl.h"
#endif

#define	CMAPSIZE	256			/* max size of colour map */
#define	WM_VOGLE	WM_USER			/* VOGLE's own window message */
#define	WM_VOGLE_ERROR	WM_USER + 1		/* VOGLE's ERROR window message */

static	HANDLE		hInst;			/* current instance */
static	HWND		hWinder;		/* handle of VOGLE's window */
static	HDC		hMemoryDC;		/* memory device context of hWinder */
static	HPEN		hPen;
static	HBITMAP		hBitmap = (HBITMAP) 0;	/* bitmap compatible with hWinder's */

static	int		CmdShow;		/* TRUE if window to be shown on start-up */
static	unsigned int	cur_color = 1;		/* index of current colour in carray */
static	COLORREF	carray[CMAPSIZE];	/* list of known colour values */
static	int		use_front_buf = TRUE;	/* FALSE if double buffering */
static	int		vogle_active = TRUE;	/* FALSE when app has been made inactive (iconized) */
static	int		have_mouse = FALSE;	/* TRUE if we have a mouse on our PC - errggh! */
static	int		xMouse = 0;		/* last known mouse position */
static	int		yMouse = 0;
static	int		LButton = 0,		/* mouse button states: 1 == down */
			MButton = 0,
			RButton = 0;
static	int		sizeX = 0;		/* width of VOGLE window */
static	int		sizeY = 0;		/* height of VOGLE window */
static	short int	xChar;			/* width of system font character */
static	short int	yChar;			/* height of system font character */
static	char		current_key = 0;	/* current key pressed at keyboard; 0 if no key is being pressed */
static	RECT		BackRect;		/* bitmap area to be swapped on next mswin_swapbuf() call */
static	int		argc;			/* command line argument count found by us */
static	char		**argv;			/* command line argument list built by us */

static	void		mswin_invalidate_rect(LPRECT);
static	void		mswin_cleanup();
static	void		swap(int *, int *);
char			*strdup(char *);

int PASCAL		WinMain(HANDLE, HANDLE, LPSTR, int);
BOOL			InitApplication(HANDLE);
BOOL			InitInstance(HANDLE, int);
long FAR PASCAL		MainWndProc(HWND, unsigned, WORD, LONG);
BOOL FAR PASCAL		About(HWND, unsigned, WORD, LONG);

/*
 * Does nothing....
 */
static
noop()
{
	return(-1);
}

/****************************************************************************

    FUNCTION: WinMain(HANDLE, HANDLE, LPSTR, int)

    PURPOSE: calls initialization function, calls user's VOGLE program (main)

    COMMENTS:

        Windows recognizes this function by name as the initial entry point 
        for the program.  This function calls the application initialization 
        routine, if no other instance of the program is running, and always 
        calls the instance initialization routine.  After breaking the command
	line passed from Windows, it executes the routine main(). This is a
	standard looking main program. It should call vinit() as soon as
	possible so that the program's window can be shown. When main
	terminates, a QUIT message is posted that starts the termination
	of the program. The program terminates in mswin_msg_loop.

        If this function must abort before entering the message loop, it 
        returns the conventional value NULL.  

****************************************************************************/

int PASCAL	WinMain(hInstance, hPrevInstance, lpCmdLine, nCmdShow)
HANDLE	hInstance;				/* current instance	     */
HANDLE	hPrevInstance;				/* previous instance	     */
LPSTR	lpCmdLine;				/* command line		     */
int	nCmdShow;				/* show-window type (open/icon) */
{
	int	i;
	LPSTR	p;
	char	v[100];

	if (!hPrevInstance)			/* Other instances of app running? */
		if (!InitApplication(hInstance))/* Initialize shared things */
			return (FALSE);		/* Exits if unable to initialize     */

	/* Perform initializations that apply to a specific instance */
	if (!InitInstance(hInstance, nCmdShow))
		return (FALSE);

	/* setup command line arguments argc, argv */
	argc = 0;
	p = lpCmdLine;
	if (p) {
		while (*p) {
			while (*p == ' ') p++;
			if (!*p) break;
			argc++;
			while ((*p != ' ') && (*p != '\0')) p++;
			}
		}

	argv = (char **)calloc(argc+1, sizeof(char *));
	if (argv) {
#ifdef VOGLE
		argv[0] = strdup("vogle");
#else
		argv[0] = strdup("vogl");
#endif
		argc = 1;
		p = lpCmdLine;
		while (*p) {
			while (*p == ' ') p++;
			if (!*p) break;
			i = 0;
			while ((*p != ' ') && (*p != '\0')) v[i++] = *(p++);
			v[i] = *p;
			argv[argc++] = strdup(v);
			}
		}
	else argc = 0;

	/* call user's program - it should call vinit() as soon as possible */
	main(argc, argv);

	/* we're done */
	PostQuitMessage(0);

	/* user's program has finished - clean up the outstanding messages */
	return (mswin_msg_loop());
}



/****************************************************************************

    FUNCTION: InitApplication(HANDLE)

    PURPOSE: Initializes window data and registers window class

    COMMENTS:

        This function is called at initialization time only if no other 
        instances of the application are running.  This function performs 
        initialization tasks that can be done once for any number of running 
        instances.  

        In this case, we initialize a window class by filling out a data 
        structure of type WNDCLASS and calling the Windows RegisterClass() 
        function.  Since all instances of this application use the same window 
        class, we only need to do this when the first instance is initialized.  


****************************************************************************/

BOOL	InitApplication(hInstance)
HANDLE	hInstance;				/* current instance	     */
{
	WNDCLASS	wc;

	/* Fill in window class structure with parameters that describe the       */
	/* main window.                                                           */

	wc.style = NULL;			/* Class style(s).                    */
	wc.lpfnWndProc = MainWndProc;		/* Function to retrieve messages for  */
						/* windows of this class.             */
	wc.cbClsExtra = 0;			/* No per-class extra data.           */
	wc.cbWndExtra = 0;			/* No per-window extra data.          */
	wc.hInstance = hInstance;		/* Application that owns the class.   */
	wc.hIcon = LoadIcon(hInstance, "VOGLEICON");
	wc.hCursor = LoadCursor(NULL, IDC_ARROW);
	wc.hbrBackground = GetStockObject(BLACK_BRUSH); 
	wc.lpszMenuName =  "VogleMenu";		/* Name of menu resource in .RC file. */
	wc.lpszClassName = "VogleWClass";	/* Name used in call to CreateWindow. */

	/* Register the window class and return success/failure code. */
	return (RegisterClass(&wc));
}



/****************************************************************************

    FUNCTION:  InitInstance(HANDLE, int)

    PURPOSE:  Saves instance handle and creates main window

    COMMENTS:

        This function is called at initialization time for every instance of 
        this application.  This function performs initialization tasks that 
        cannot be shared by multiple instances.  

        In this case, we save the instance handle in a static variable and 
        create and display the main program window.  
        
****************************************************************************/

BOOL	InitInstance(hInstance, nCmdShow)
HANDLE	hInstance;				/* Current instance identifier.       */
int	nCmdShow;				/* Param for first ShowWindow() call. */
{
	/* Save the instance handle in static variable, which will be used in  */
	/* many subsequent calls from this application to Windows.            */
	/* Also save the handle of the main window - hWinder */

	hInst = hInstance;
	CmdShow = nCmdShow;

	return (TRUE);
}


/****************************************************************************

    FUNCTION: MainWndProc(HWND, unsigned, WORD, LONG)

    PURPOSE:  Processes messages

    MESSAGES:

	WM_PAINT	- copy bit image from back buffer to front buffer
	WM_MOUSEMOVE	- track position of mouse as it moves over window
	WM_LBUTTONDOWN	- track state of mouse button
	WM_LBUTTONUP	- track state of mouse button
	WM_MBUTTONDOWN	- track state of mouse button
	WM_MBUTTONUP	- track state of mouse button
	WM_RBUTTONDOWN	- track state of mouse button
	WM_RBUTTONUP	- track state of mouse button
	WM_CHAR		- track last key pressed on keyboard
	WM_KEYUP	- track last key pressed on keyboard
	WM_ACTIVATE	- track active state of this application
	WM_COMMAND	- application menu (About dialog box)
	WM_DESTROY	- destroy window

    COMMENTS:

	To process the IDM_ABOUT message, call MakeProcInstance() to get the
	current instance address of the About() function.  Then call Dialog
	box which will create the box according to the information in your
	Vogle.rc file and turn control over to the About() function.	When
	it returns, free the intance address.

****************************************************************************/

long FAR PASCAL MainWndProc(hWnd, message, wParam, lParam)
HWND	hWnd;					/* window handle		     */
unsigned message;				/* type of message		     */
WORD	wParam;					/* additional information	     */
LONG	lParam;					/* additional information	     */
{
	HDC		hDC;
	HBITMAP		hOldBitmap;
	PAINTSTRUCT	ps;

	FARPROC lpProcAbout;			/* pointer to the "About" function */

	switch (message) {

		case WM_PAINT :
			if (!hBitmap) {
				/* no back buffer yet - let's use default procedure */
				return (DefWindowProc(hWnd, message, wParam, lParam));
				}
			BeginPaint(hWnd, &ps);
			hOldBitmap = SelectObject(hMemoryDC, hBitmap);
			hDC = GetDC(hWinder);
			BitBlt(hDC,
				ps.rcPaint.left, ps.rcPaint.top,
				ps.rcPaint.right - ps.rcPaint.left,
				ps.rcPaint.bottom - ps.rcPaint.top,
				hMemoryDC,
				ps.rcPaint.left, ps.rcPaint.top,
				SRCCOPY);
			SelectObject(hMemoryDC, hOldBitmap);
			ReleaseDC(hWnd, hDC);
			EndPaint(hWnd, &ps);
			break;

		case WM_MOUSEMOVE :
			xMouse = LOWORD(lParam);
			yMouse = HIWORD(lParam);
			break;

		case WM_LBUTTONDOWN :	LButton = 1; break;
		case WM_LBUTTONUP :	LButton = 0; break;
		case WM_MBUTTONDOWN :	MButton = 1; break;
		case WM_MBUTTONUP :	MButton = 0; break;
		case WM_RBUTTONDOWN :	RButton = 1; break;
		case WM_RBUTTONUP :	RButton = 0; break;

		case WM_CHAR:
			/* user has pressed a character - save it for mswin_checkkey() */
			current_key = wParam;
			break;

		case WM_KEYUP:
			/* user has released a key - if it is a meta-key (SHIFT, CTRL) ignore it */
			/* otherwise we assume that the key released was the one that generated  */
			/* the WM_CHAR message and we reset current_key. This assumption is      */
			/* probably VERY wrong!                                                  */
			switch (wParam) {
				case VK_SHIFT :
				case VK_CONTROL :
				case VK_MENU :
				case VK_CAPITAL :
				case VK_NUMLOCK :
					break;
				default : current_key = 0;
				}
			break;

		case WM_ACTIVATE :
			vogle_active = !IsIconic(hWnd);
			break;

		case WM_VOGLE_ERROR :
			/* VOGLE error occurred - display error message in a window, then QUIT */
			MessageBox(hWinder, (LPSTR)lParam, "A VOGLE error has occurred!", MB_APPLMODAL | MB_ICONEXCLAMATION |  MB_OK);
			PostMessage(hWinder, WM_DESTROY, NULL, NULL);
			break;

		case WM_COMMAND :	   /* message: command from application menu */
			if (wParam == IDM_ABOUT) {
				lpProcAbout = MakeProcInstance(About, hInst);
				DialogBox(hInst,	/* current instance	     */
					"AboutBox",	/* resource to use	     */
					hWnd,		/* parent handle	     */
					lpProcAbout);	/* About() instance address */
				FreeProcInstance(lpProcAbout);
				}
			else return (DefWindowProc(hWnd, message, wParam, lParam));
			break;

		case WM_DESTROY :			/* message: window being destroyed */
		    	PostQuitMessage(0);
		    	break;

		default :				/* Passes it on if unproccessed    */
	    		return (DefWindowProc(hWnd, message, wParam, lParam));
		}
	return (NULL);
}


/****************************************************************************

    FUNCTION: About(HWND, unsigned, WORD, LONG)

    PURPOSE:  Processes messages for "About" dialog box

    MESSAGES:

	WM_INITDIALOG - initialize dialog box
	WM_COMMAND    - Input received

    COMMENTS:

	No initialization is needed for this particular dialog box, but TRUE
	must be returned to Windows.

	Wait for user to click on "Ok" button, then close the dialog box.

****************************************************************************/

BOOL FAR PASCAL About(hDlg, message, wParam, lParam)
HWND	hDlg;					/* window handle of the dialog box */
unsigned message;				/* type of message                 */
WORD	wParam;					/* message-specific information    */
LONG	lParam;
{
	switch (message) {
		case WM_INITDIALOG:		/* message: initialize dialog box */
			return (TRUE);

		case WM_COMMAND:		/* message: received a command */
			if (wParam == IDOK ||	/* "OK" box selected?	     */
			    wParam == IDCANCEL) {	/* System menu close command? */
				EndDialog(hDlg, TRUE);	/* Exits the dialog box	     */
				return (TRUE);
				}
		break;
		}
	return (FALSE);				/* Didn't process a message    */
}



/****************************************************************************
	VOGLE bits
****************************************************************************/

static
int		mswin_msg_loop()
{
	MSG	msg;

	/* put our own message in the message queue */
	PostMessage(hWinder, WM_VOGLE, NULL, NULL);

	while (GetMessage(&msg, NULL, NULL, NULL) || !vogle_active) {
		TranslateMessage(&msg);	   /* Translates virtual key codes	     */
		DispatchMessage(&msg);	   /* Dispatches message to window	     */
		if (msg.message == WM_VOGLE) break;
		}

	if (msg.message == WM_QUIT) {
		mswin_cleanup();
		exit(msg.wParam);
		}
	return (msg.wParam);
}



static
void		mswin_cleanup()
{
	int	i;

	DeleteObject(hBitmap);
	DeleteDC(hMemoryDC);
	for (i = 0; i < argc; i++) free(argv[i]);
	free(argv);
}



static
void		mswin_invalidate_rect(r)
	LPRECT	r;
{
	if (r->left < BackRect.left) BackRect.left = r->left;
	if (r->right > BackRect.right) BackRect.right = r->right;
	if (r->top < BackRect.top) BackRect.top = r->top;
	if (r->bottom > BackRect.bottom) BackRect.bottom = r->bottom;
}



int		mswin_init()
{
	RECT	rect;
	HDC	hDC;				/* device context of Desktop Window */
	TEXTMETRIC	tm;
	RECT	r;
	HBRUSH	hBrush,
		hOldBrush;
	HBITMAP	hOldBitmap;
	int	prefx, prefy,			/* values for the preferred */
		prefxs, prefys,			/*   window size and position */
		x, y;

	char	buf[128];


	/* determine a preferred size and position */
	getprefposandsize(&prefx, &prefy, &prefxs, &prefys);
	if (prefx > -1) {
		x = prefx;
		y = prefy;
		}
	else {
		x = y = CW_USEDEFAULT;
		}

	if (prefxs > -1) {
		sizeX = prefxs;
		sizeY = prefys;
		}
	else {
		/* figure out a good size based on display device capabilities */

		hDC = GetDC(GetDesktopWindow());
		sizeX = sizeY = ((GetDeviceCaps(hDC, HORZRES) * 2) / 3);
		ReleaseDC(GetDesktopWindow(), hDC);
		}


	/* Create a main window for this application instance */
#ifdef VOGLE
	strcpy(buf, "Vogle Application");
#else
	if (vdevice.wintitle)
		strcpy(buf, vdevice.wintitle);
	else
		strcpy(buf, "VOGL Application");
#endif

	hWinder = CreateWindow(
#ifdef VOGLE
		"VogleWClass",			/* See RegisterClass() call.          */

#else
		"VogleWClass",			/* See RegisterClass() call.          */
#endif
		buf,				/* Text for window title bar.         */
		WS_OVERLAPPEDWINDOW,            /* Window style.                      */
		x,                  		/* Default horizontal position.       */
		y,				/* Default vertical position.         */
		sizeX,				/* Default width.                     */
		45 + sizeY,			/* Default height.  + a bit                   */
		NULL,                           /* Overlapped windows have no parent. */
		NULL,                           /* Use the window class menu.         */
		hInst,                          /* This instance owns this window.    */
		NULL                            /* Pointer not needed.                */
		);

	/* If window could not be created, return "failure" */
	if (!hWinder) return (FALSE);
	hDC = GetDC(hWinder);
	SetViewportOrg(hDC, 0, 0);
	ShowWindow(hWinder, CmdShow);	/* Show the window                        */
	UpdateWindow(hWinder);          /* Sends WM_PAINT message                 */

	mswin_msg_loop();
 
	/*
	 *  Let VOGLE know about the dimensions of the window known as hWinder.
	 */

	hDC = GetDC(GetDesktopWindow());
	vdevice.depth = GetDeviceCaps(hDC, PLANES);
	ReleaseDC(GetDesktopWindow(), hDC);

	/* determine size of default (system) font */
	hDC = GetDC(hWinder);
	GetTextMetrics(hDC, &tm);
	xChar = tm.tmMaxCharWidth;
	yChar = tm.tmHeight + tm.tmExternalLeading;

	/* set window size - make it square */
	vdevice.sizeX = vdevice.sizeY = min(sizeX, sizeY) - 1;
 	vdevice.sizeSx = sizeX - 1;
	vdevice.sizeSy = sizeY - 1;

	/* create a bitmap and memory context compatible with hWinder */
	hMemoryDC = CreateCompatibleDC(hDC);
	hBitmap = CreateCompatibleBitmap(hDC, vdevice.sizeSx, vdevice.sizeSy);

	/* clear the bitmap */
	hOldBitmap = SelectObject(hMemoryDC, hBitmap);
	hBrush = GetStockObject(BLACK_BRUSH);
	if (hBrush) {
		hOldBrush = SelectObject(hMemoryDC, hBrush);
		r.top = r.left = 0;
		r.right = vdevice.sizeSx;
		r.bottom = vdevice.sizeSy;
		FillRect(hMemoryDC, &r, hBrush);
		SelectObject(hMemoryDC, hOldBrush);
		}
	SelectObject(hMemoryDC, hOldBitmap);

	ReleaseDC(hWinder, hDC);

	/* set up colour indices */
	carray[0] = RGB(  0,   0,   0);
	carray[1] = RGB(255,   0,   0);
	carray[2] = RGB(  0, 255,   0);
	carray[3] = RGB(255, 255,   0);
	carray[4] = RGB(  0,   0, 255);
	carray[5] = RGB(255,   0, 255);
	carray[6] = RGB(  0, 255, 255);
	carray[7] = RGB(255, 255, 255);

	have_mouse = GetSystemMetrics(SM_MOUSEPRESENT);

	mswin_msg_loop();


	return TRUE;
}



int		mswin_exit()
{
	PostQuitMessage(0);
	mswin_msg_loop();
	return 1;
}



int		mswin_clear()
{
	HDC	hDC;
	RECT	r;
	HBRUSH	hBrush,
		hOldBrush;
	HBITMAP	hOldBitmap;

	hOldBitmap = SelectObject(hMemoryDC, hBitmap);
	GetClientRect(hWinder, &r);
	hBrush = CreateSolidBrush(carray[cur_color]);
	if (hBrush) {
		hOldBrush = SelectObject(hMemoryDC, hBrush);
		FillRect(hMemoryDC, &r, hBrush);
		if (use_front_buf) {
			hDC = GetDC(hWinder);
			BitBlt(hDC,
				r.left, r.top,
				r.right - r.left,
				r.bottom - r.top,
				hMemoryDC,
				r.left, r.top,
				SRCCOPY);
			ReleaseDC(hWinder, hDC);
			}
		else {
			/* invalidate the region in front buffer */
			mswin_invalidate_rect(&r);
			}
		SelectObject(hMemoryDC, hOldBrush);
		DeleteObject(hBrush);
		}
	SelectObject(hMemoryDC, hOldBitmap);

	mswin_msg_loop();
}



int		mswin_draw(x, y)
int		x, y;
{
	HDC	hDC;
	HPEN	hOldPen;
	HBITMAP	hOldBitmap;
	RECT	r;

	hOldBitmap = SelectObject(hMemoryDC, hBitmap);
	hPen = CreatePen(PS_DOT, 3, carray[cur_color]);
	if (hPen) {
		hOldPen = SelectObject(hMemoryDC, hPen);
		MoveTo(hMemoryDC, r.left = vdevice.cpVx, r.top = vdevice.sizeSy - vdevice.cpVy);
		LineTo(hMemoryDC, r.right = x, r.bottom = vdevice.sizeSy - y);
		if (use_front_buf) {
			if (r.top > r.bottom) swap(&(r.top), &(r.bottom));
			if (r.left > r.right) swap(&(r.left), &(r.right));
			hDC = GetDC(hWinder);
			BitBlt(hDC,
				r.left, r.top,
				r.right - r.left + 1,
				r.bottom - r.top + 1,
				hMemoryDC,
				r.left, r.top,
				SRCCOPY);
			ReleaseDC(hWinder, hDC);
			}
		else {
			/* invalidate the region in front buffer */
			mswin_invalidate_rect(&r);
			}
		SelectObject(hMemoryDC, hOldPen);
		DeleteObject(hPen);
		}
	hOldBitmap = SelectObject(hMemoryDC, hOldBitmap);

	mswin_msg_loop();
}



int		mswin_fill(n, x, y)
int		n;
int		*x, *y;
{
	HDC	hDC;
	LPPOINT	lpPoint;
	int	OldPolyMode;
	int	i;
	HBRUSH	hBrush,
		hOldBrush;
	HBITMAP	hOldBitmap;
	RECT	r;

	hOldBitmap = SelectObject(hMemoryDC, hBitmap);

	if (lpPoint = (LPPOINT) malloc(sizeof(POINT) * n)) {
 		hBrush = CreateSolidBrush(carray[cur_color]);
		if (hBrush) {
			r.left = r.right = x[0];
			r.top = r.bottom = y[0];
			for (i = 0; i < n; i++) {
				lpPoint[i].x = x[i];
				if (lpPoint[i].x < r.left) r.left = lpPoint[i].x;
				else if (lpPoint[i].x > r.right) r.right = lpPoint[i].x;

				lpPoint[i].y = vdevice.sizeSy - y[i];
				if (lpPoint[i].y < r.top) r.top = lpPoint[i].y;
				else if (lpPoint[i].y > r.bottom) r.bottom = lpPoint[i].y;
				}
			hOldBrush = SelectObject(hMemoryDC, hBrush);
			OldPolyMode = GetPolyFillMode(hMemoryDC);
			SetPolyFillMode(hMemoryDC, WINDING);
			Polygon(hMemoryDC, lpPoint, n);
			if (use_front_buf) {
				hDC = GetDC(hWinder);
				BitBlt(hDC,
					r.left, r.top,
					r.right - r.left + 1,
					r.bottom - r.top + 1,
					hMemoryDC,
					r.left, r.top,
					SRCCOPY);
				ReleaseDC(hWinder, hDC);
				}
			else {
				/* invalidate the region in front buffer */
				mswin_invalidate_rect(&r);
				}

			SetPolyFillMode(hMemoryDC, OldPolyMode);
			SelectObject(hMemoryDC, hOldBrush);
			DeleteObject(hBrush);
			}
		free(lpPoint);
		}
	hOldBitmap = SelectObject(hMemoryDC, hOldBitmap);

	vdevice.cpVx = x[n-1];
	vdevice.cpVy = y[n-1];

	mswin_msg_loop();
}



int		mswin_color(i)
int		i;
{
	cur_color = (unsigned)i;
}



int		mswin_mapcolor(i, r, g, b)
int		i, r, g, b;
{
	if (i >= CMAPSIZE) return(-1);

	carray[i] = RGB(r, g, b);

	mswin_msg_loop();
	return 0;
}



int		mswin_char(c)
char	c;
{
	HDC	hDC;		
	char	s[2];
	HBITMAP	hOldBitmap;
	RECT	r;
	int	OldBkMode;
	DWORD	OldTextColour;

	hOldBitmap = SelectObject(hMemoryDC, hBitmap);

	s[0] = c;
	s[1] = '\0';

	OldBkMode = SetBkMode(hMemoryDC, TRANSPARENT);
	OldTextColour = SetTextColor(hMemoryDC, carray[cur_color]);
	TextOut(hMemoryDC, r.left = vdevice.cpVx, r.top = vdevice.sizeSy - vdevice.cpVy, s, 1);
	SetTextColor(hMemoryDC, OldTextColour);
	SetBkMode(hMemoryDC, OldBkMode);

	if (use_front_buf) {
		hDC = GetDC(hWinder);
		BitBlt(hDC,
			r.left, r.top,
			xChar,
			yChar,
			hMemoryDC,
			r.left, r.top,
			SRCCOPY);
		ReleaseDC(hWinder, hDC);
		}
	else {
		/* invalidate the region in front buffer */
		r.right = r.left + xChar;
		r.bottom = r.top + yChar;
		mswin_invalidate_rect(&r);
		}

	hOldBitmap = SelectObject(hMemoryDC, hOldBitmap);

	mswin_msg_loop();
}



int		mswin_string(s)
char	*s;
{
	HDC	hDC;
	HBITMAP	hOldBitmap;
	RECT	r;
	int	n;
	int	OldBkMode;
	DWORD	OldTextColour;

	hOldBitmap = SelectObject(hMemoryDC, hBitmap);

	OldBkMode = SetBkMode(hMemoryDC, TRANSPARENT);
	OldTextColour = SetTextColor(hMemoryDC, carray[cur_color]);
	TextOut(hMemoryDC, r.left = vdevice.cpVx, r.top = vdevice.sizeSy - vdevice.cpVy, s, n = strlen(s));
	SetTextColor(hMemoryDC, OldTextColour);
	SetBkMode(hMemoryDC, OldBkMode);
	if (use_front_buf) {
		hDC = GetDC(hWinder);
		BitBlt(hDC,
			r.left, r.top,
			n * xChar,
			yChar,
			hMemoryDC,
			r.left, r.top,
			SRCCOPY);
		ReleaseDC(hWinder, hDC);
		}
	else {
		/* invalidate the region in front buffer */
		r.right = r.left + n * xChar;
		r.bottom = r.top + yChar;
		mswin_invalidate_rect(&r);
		}

	hOldBitmap = SelectObject(hMemoryDC, hOldBitmap);

	mswin_msg_loop();
}



int		mswin_font(font)
char	*font;
{
	/* NOT IMPLEMENTED */
	mswin_msg_loop();
	return 1;
}



int		mswin_frontbuf()
{
	use_front_buf = TRUE;
	mswin_msg_loop();
}



int		mswin_backbuf()
{
	use_front_buf = FALSE;
	BackRect.left = BackRect.top = 32765;
	BackRect.right = BackRect.bottom = -32765;
	mswin_msg_loop();
	return 1;
}



int		mswin_swapbuf()
{
	HDC	hDC;
	RECT	r;
	HBITMAP	hOldBitmap;

	hDC = GetDC(hWinder);
	hOldBitmap = SelectObject(hMemoryDC, hBitmap);
#ifdef GUNGE
	GetClientRect(hWinder, &r);
	BitBlt(hDC,
		r.left, r.top,
		r.right - r.left,
		r.bottom - r.top,
		hMemoryDC,
		r.left, r.top,
		SRCCOPY);
#endif
	BitBlt(hDC,
		BackRect.left, BackRect.top,
		BackRect.right - BackRect.left,
		BackRect.bottom - BackRect.top,
		hMemoryDC,
		BackRect.left, BackRect.top,
		SRCCOPY);
	BackRect.left = BackRect.top = 32765;
	BackRect.right = BackRect.bottom = -32765;

	ReleaseDC(hWinder, hDC);
	SelectObject(hMemoryDC, hOldBitmap);

	mswin_msg_loop();
}



int		mswin_checkkey()
{
	MSG	Msg;

	while (PeekMessage(&Msg, NULL, NULL, NULL, PM_NOREMOVE)) {
		GetMessage(&Msg, NULL, NULL, NULL);
		TranslateMessage(&Msg);	   /* Translates virtual key codes	     */
		DispatchMessage(&Msg);	   /* Dispatches message to window	     */
		}
	if (Msg.message == WM_QUIT) {
		mswin_cleanup();
		exit(Msg.wParam);
		}

	return current_key;
}



int		mswin_getkey()
{
	MSG	Msg;
	WORD	quit;

	while (GetMessage(&Msg, NULL, NULL, NULL)) {
		TranslateMessage(&Msg);	   /* Translates virtual key codes	     */
		DispatchMessage(&Msg);	   /* Dispatches message to window	     */
		if (Msg.message == WM_CHAR && vogle_active) break;
		}
	if (Msg.message == WM_QUIT) {
		mswin_cleanup();
		exit(Msg.wParam);
		}
	return (Msg.wParam);
}



int		mswin_locator(x, y)
int		*x, *y;
{
	mswin_msg_loop();

	*x = xMouse;
	*y = vdevice.sizeSy - yMouse;
	return (LButton | (MButton << 2) | (RButton << 1));
}



int		mswin_verror(msg)
char		*msg;
{
	PostMessage(hWinder, WM_VOGLE_ERROR, NULL, (DWORD)((LPSTR) msg));
	mswin_msg_loop();
}


#ifndef VOGLE
/*
 * Haven't got around to these for VOGL yet.
 */
void
mswin_setls()
{}

void
mswin_setlw()
{}
#endif

/*
 * the device entry
 */
static DevEntry mswindev = {
	"mswin",
	"large",
	"small",
	mswin_backbuf,
	mswin_char,
	mswin_checkkey,
	mswin_clear,
	mswin_color,
	mswin_draw,
	mswin_exit,
	mswin_fill,
	mswin_font,
	mswin_frontbuf,
	mswin_getkey,
	mswin_init,
	mswin_locator,
	mswin_mapcolor,
#ifndef VOGLE
	mswin_setls,
	mswin_setlw,
#endif
	mswin_string,
	mswin_swapbuf,
	noop
};

/*
 * _mswin_devcpy
 *
 *	copy the mswin device into vdevice.dev.
 */
int		_mswin_devcpy()
{
	vdevice.dev = mswindev;
}



static
void	swap(a, b)
int		*a;
int		*b;
{
	int	t;

	t = *a;
	*a = *b;
	*b = t;
}
