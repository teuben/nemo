#ifndef _Xtc_h
#define _Xtc_h

/*
 *  Xtc header
 */

/*#define SCREEN_WIDTH 640
#define SCREEN_HEIGHT 400
*/
#define SCREEN_WIDTH 512
#define SCREEN_HEIGHT 512


enum color_macros {
	BLACK, BLUE, GREEN, CYAN, RED, MAGENTA, BROWN, LIGHTGRAY,
	DARKGRAY, LIGHTBLUE, LIGHTGREEN, LIGHTCYAN, LIGHTRED, LIGHTMAGENTA,
	YELLOW, WHITE, OLD_NCOLORS, };

#define NCOLORS 32

enum graphics_functions {
	COPY_PUT, XOR_PUT, OR_PUT, AND_PUT, NOT_PUT, };

enum line_styles {
	SOLID_LINE, DOTTED_LINE, CENTER_LINE, DASHED_LINE, };

enum line_widths { NORM_WIDTH = 1, THICK_WIDTH = 3, };

void arc(), bar(), circle(), cleardevice(), closegraph(), drawpoly(),
     ellipse(), fillellipse(), fillpoly(), getimage(), 
     initgraph(), line(), linerel(), lineto(), moverel(), moveto(),
     outtextxy(), pieslice(), putimage(), putpixel(), rectangle(),
	 sector(), setbkcolor(), setlinestyle(), setrgbpalette(),
	 setwritemode();

int getbkcolor(), getcolor(), getmaxcolor(), getmaxx(), getmaxy(),
	getpixel(), getx(), gety(), imagesize(),
	textheight(), textwidth();
    
int cur_pointer_x, cur_pointer_y;

void xtcmainloop(), xflush();
int	xgetbutton();

void copy_from_pixmap(void);

void setcolor_pixmap(),fillellipse_pixmap(),cleardevice_pixmap();
void circle_pixmap(),outtextxy_pixmap();
void rectangle_pixmap();

#endif /* _Xtc_h */
