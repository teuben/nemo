/*
 *   MOVIE_VS:  show a number of pixrect (raster) files as a movie.
 *
 *		This version is specifically designed to work under
 *		suntools in the sunview environment, and can also
 *		be compiler for the NEMO environment (-DNEMO)
 *  
 *      This version of MOVIE is using command mode thorugh the keyboard
 *
 *      q, ^C       quit
 *      h, ?       help menu
 *      ssf, +      single step forward
 *      ssb, -      single step backward
 *      ff          fast forward
 *      fb          fast backward
 *
 *	NOTE: a few names of routines had to be changed because sunview
 *		header files do not allow certain names
 *
 *
 *
 *		dec-88	created from movie -- Peter Teuben
 *	     26-sep-89  bug removed for SUN OS 4.0 (sun 3/4)
 *           27-nov-89  fix select() upgrade for SUN OS 4.0
 *           12-mar-91  working on adding -DNEMO nemo_main() section 
 *			removed coredump section when frame.0 not present
 */

#define WINDOWS 1
#define SUNOS4  1

#include <stdio.h>
#if defined(WINDOWS)
#  include <suntool/sunview.h>
#  include <suntool/canvas.h>
#else
#  include <pixrect/pixrect_hs.h>
#endif

#include <stdio.h>
#include <sgtty.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/time.h>

#define NULLPR (struct pixrect *)NULL


Frame base_frame;
Canvas canvas;
Pixwin *pw;

struct pixrect *screen;
struct pixfont *pf_default(), *pfont;
int xsize, ysize;

int oldr, oldc, row, col, oldx, oldy, sx, sy, ix, iy, go;
int oldm, mrow, oldmy, my, imy, page, px, oldpx, pc, oldpc;
int bl, bc, br;
int MouseFD;
struct pr_prpos    xywh;
int (*Button[4])();
int (*MenuFunction[16][32])();

FILE *mouse;

#define LEFT 1
#define CENTRE 2
#define RIGHT 4
#define L 1
#define C 2
#define R 3
#define MenuX1 910
#define MenuY1 60
#define MenuY2 764
#define MenuWidth 240
#define MenuHeight 22
#define MenuPages 16
#define MenuRows 32
#define LastMenuY MenuY2 - MenuHeight

#define PageP1 MenuX1
#define PageWidth 30
#define PageP2 PageP1 + PageWidth
#define PageEnd PageP2 + PageWidth

#define MLineX 910
#define MLineY 30
#define MLWidth 80

#define PIX_XOR (PIX_SRC|PIX_DST)&PIX_NOT(PIX_SRC&PIX_DST)
#define MAXFRAMES 256
struct pixrect *movie[MAXFRAMES];
int nframe, lastframe, frame, rolling, screenframe, framerate;
struct timeval delay;
int delaytab_s[201], delaytab_us[201];
int framex, framey, frameop;

char menu[MenuPages][MenuRows+1][128], MouseString[4][11];
int MainX1, MainY1, MainX2, MainY2;

char ftext[MAXFRAMES][80];  /* label for frame */
char scriptfile[128], basename[128], firstfile[128];

#if !defined(NEMO)

main(argc,argv)
int argc;
char *argv[];
{
    int argn;

    scriptfile[0] == '\0';
    strcpy(basename,"frame");

    for (argn = 1; argn < argc; argn++) {
        if (argv[argn][0] == '-') {
            switch (argv[argn][1]) {
                case 'p': /* frame position */
                    framey = atoi(argv[++argn]);
                    framex = atoi(argv[++argn]);
                    break;
                case 's': /* script */
                    strcpy(scriptfile,argv[++argn]);
                    break;
		case 'f': /* basename */
	            strcpy(basename,argv[++argn]);
                    break;
                default:
                    Usage();
                    exit(0);
            }
        } else {
            Usage();
            exit(0);
        }
    }

    if (Movie_Init())       /* if things to do */
        Movie_Menu();       /* do them */
    Movie_Fini();           /* clear up any mess done */
    exit(0);
}
#else
#include <stdinc.h>
#include <getparam.h>

string defv[] = {
    "frame=frame\n      Basename of frame files",
    "file=\n            File extensions to use [0,1,...]",
    "script=\n          Or: use a scriptfile",
    "VERSION=1.0\n      12-mar-91 PJT",
    NULL,
};

nemo_main()
{
    string frame, script;
    int i, nframe, iframe[MAXFRAMES];
    stream fp;
    bool   delete=FALSE;

    script = getparam("script");
    if (*script) {              /* use script file */
        strcpy(scriptfile,script);
    } else {
        frame = getparam("frame");
        nframe = nemoinpi(getparam("file"), iframe, MAXFRAMES);
        if (nframe==0) {
            strcpy(basename,frame);
        } else {
            delete = TRUE;      /* use temporary file */
            sprintf(scriptfile,"/tmp/movie_sv%d",getpid());
            fp = stropen(scriptfile,"w!");      /* ALWAYS write it */
            for (i=0; i<nframe; i++)
                fprintf(fp,"%s.%d\n",frame,iframe[i]);
            strclose(fp);
        }
    }

    if (Movie_Init())       /* if things to do */
        Movie_Menu();       /* do them */
    Movie_Fini();           /* clear up any mess done */

    if (delete) unlink(scriptfile);     /* delete temp script file */
}
#endif

Usage()
{
    fprintf(stderr,"movie [-p %%d %%d] [-s %%s] [-f %%s]\n");
    fprintf(stderr,"    -p X Y  Frame row col position\n");
    fprintf(stderr,"    -s S    Script file name\n");
    fprintf(stderr,"    -f B    Basename for frame-files \n");
}

Movie_Menu()
{
    int dx, dy, but;
#if 0
    DrawMenuCursor(px,my);

    while (go) {
        ReadMouse(&dx,&dy,&but);
        MoveMenuCursor(dx,dy);
        if (but) MenuHandle(my,but);
    }

    DrawMenuCursor(px,my);
#else
    printf ("No menu yet\n");
#endif
}

MoveMenuCursor(dx,dy)
int dx,dy;
{
    oldm = mrow;
    oldmy = my;
    oldpc = pc;
    oldpx = px;
    imy -= dy;
    px += dx;
    if (imy < MenuY1) imy = MenuY1;
    if (imy > MenuY2) imy = MenuY2;
    if (px < PageP1) px = PageP1;
    if (px > PageEnd) px = PageEnd;

    mrow = 1 + (imy - MenuY1) / MenuHeight;
    my = MenuY1 + (mrow - 1) * MenuHeight;
    pc = (px > PageP2);

    if (oldm != mrow || oldpc != pc) {
        DrawMenuCursor(oldpx,oldmy);
	RefreshMenu();
        DrawMenuCursor(px,my);
    }
}

DrawMenuCursor(x,y)
int x,y;
{
    int r;

    if (screen==NULL) return;
    r = (y - MenuY1) / MenuHeight;

    if (y > LastMenuY) {
        if (x > PageP2) pr_rop(screen,PageP2,y,PageWidth,MenuHeight,
                PIX_NOT(PIX_DST),NULLPR,0,0);
        else pr_rop(screen,PageP1,y,PageWidth,MenuHeight,
                PIX_NOT(PIX_DST),NULLPR,0,0);
    }
    else {
        pr_rop(screen,MenuX1,y,MenuWidth-1,MenuHeight,
            PIX_NOT(PIX_DST),NULLPR,0,0);

        if (menu[page][r][0] == '\0')
            pr_rop(screen,MenuX1+3,y+3,MenuWidth-7,MenuHeight-6,
                PIX_NOT(PIX_DST),NULLPR,0,0);
    }
}

BindButton(b,str,f)
int b;
char *str;
int (*f)();
{
    Button[b] = f;
    MouseLine(b,str);
}

BOGUS()
{
}

MouseLine(b,str)
int b;
char *str;
{
    MouseString[b][0] = ' ';
    strncpy(MouseString[b]+1,str,8);
    MouseString[b][9] = ' ';
    MouseString[b][10] = '\0';
    DrawMouseItem(b);
}

RefreshMouseLine()
{
    struct pr_prpos wh;

    wh.pr = screen;
    wh.pos.x = MLineX;
    wh.pos.y = MLineY + MenuHeight - 3;
    pr_rop(screen,wh.pos.x,wh.pos.y-14,MenuWidth+1,19,PIX_CLR,NULLPR,0,0);
    wh.pos.x += 8;
    wh.pos.y -= 3;
    pf_text(wh,PIX_NOT(PIX_SRC),pfont," left   ");
    wh.pos.x += MLWidth;
    pf_text(wh,PIX_NOT(PIX_SRC),pfont," centre ");
    wh.pos.x += MLWidth;
    pf_text(wh,PIX_NOT(PIX_SRC),pfont," right  ");

    DrawMouseItem(L);
    DrawMouseItem(C);
    DrawMouseItem(R);
}

DrawMouseItem(b)
int b;
{
    struct pr_prpos wh;

    wh.pr = screen;
    wh.pos.x = MLineX + (b - 1) * MLWidth;
    wh.pos.y = MLineY;
    pr_rop(screen,wh.pos.x,wh.pos.y-14,MLWidth+1,19,PIX_CLR,NULLPR,0,0);
    pf_text(wh,PIX_SRC,pfont,MouseString[b]);
}

Movie_Init()
{
    FILE *fopen();
    int one = 1;
    int i, s, us;
    double d, dd, ds, dus;
    struct rasterfile rh;
    FILE *fp, *fopen();

    if (LoadFrames() == 0) return(0);


#if defined(WINDOWS)
    printf ("WINDOWS-version:\n");      /* Check using First Frame */
    if ((fp=fopen(firstfile,"r"))==NULL) {
        printf("Error opening first frame %s\n",firstfile);
        return(0);
    }
    if (pr_load_header(fp,&rh)) {
        printf("Error in leading first rasterfile header\n");
        return(0);
    } else {
        fclose(fp);
        printf ("First rasterfile %s has npx=%d npy=%d npz=%d\n",
                firstfile, rh.ras_width, rh.ras_height, rh.ras_depth);
	xsize = rh.ras_width;
	ysize = rh.ras_height;
    }
#else
    screen = pr_open("/dev/fb");
    mouse = fopen("/dev/mouse","r");
    if (mouse == NULL) {
        (void)printf("can't open /dev/mouse\n");
        return(0);
    }
    MouseFD = fileno(mouse);
#endif

    go = 1;
    framerate = 90;
    screenframe = -1;
    frame = 0;
    framex = 0;
    framey = 0;
    frameop = PIX_SRC;

    delaytab_s [  0] = delaytab_s [200] = 0;
    delaytab_us[  0] = delaytab_us[200] = 0;
    d = 0.75;
    dd = d  / 99.0;

    for (i = 1; i < 100; i++) {
        s = (int)d;
        dus = (d - (double)s) * 100000.0;
        us = (int)dus;
        delaytab_s[i + 100] = delaytab_s[100 - i] = s;
        delaytab_us[i + 100] = delaytab_us[100 - i] = 10*us;
        d -= dd;
#if 0
        printf ("%4d => %d %d     %4d => %d %d\n",
        i+100,  delaytab_s[i + 100], delaytab_us[i + 100], 
        100-i,  delaytab_s[100 - i], delaytab_us[100 - i]);
#endif
    }

    if (1) {
        printf ("BOGUS movie\n");
        Bogus_Movie();
        return(0);
    }

    MainX1 = 0;
    MainY1 = 0;
    MainX2 = 899;
    MainY2 = 899;
    imy = MenuY1;
    my = MenuY1;
    mrow = 1;
    px = PageP1;
    pc = 0;

    xywh.pr = screen;
    xywh.pos.x = 910;
    xywh.pos.y = 10;

    page = 1;

    pfont = pf_default();
    EraseScreen();
    InitMenu();
    page = 0;
    DisplayFrame();
    RefreshMenu();
    Button[0] = BOGUS;
    BindButton(L,"Select",BOGUS);
    BindButton(C," ",BOGUS);
    BindButton(R," ",BOGUS);
    RefreshMouseLine();
    return(1);
}

LoadFrames()
{
    int ifile;
    char str[256], *malloc(), *cmap;
    FILE *fp, *fopen(), *script;
    struct pixrect *pr_load();
    int first=1;
/*    struct colormap_t *cmap;	*/		/* KLUDGE AWAY */

    firstfile[0] = '\0';            /* clear name of first frame file */
    fprintf(stderr,"Loading frames with basename %s\n",basename);
    cmap = malloc(5*4);				/* PJT KLUDGE for Sun OS 4 */
    if (scriptfile[0] == '\0') {
        for (frame=0, ifile=1; frame<MAXFRAMES && ifile; frame++) {
            sprintf(str,"%s.%d",basename,frame);
            fp = fopen(str,"r");
            if (fp == NULL) 
                ifile = 0;
            else {
                if (first) {
                    strcpy(firstfile,str);
                    first = 0;
                }
                movie[frame] = pr_load(fp,cmap);
                fclose(fp);
                fprintf(stderr,"%2d ",frame);
                if (9 == (frame % 10)) fprintf(stderr,"\n");
                fflush(stderr);
            }
        }
    } else {
        script = fopen(scriptfile,"r");
        if (script == NULL) {
            fprintf(stderr,"movie: can't read script file\"%s\"\n",scriptfile);
            return(0);
        }
        frame = 0;
        while (fscanf(script,"%s",str) != EOF) {
            fp = fopen(str,"r");
            if (fp == NULL) {
                fprintf(stderr,"movie: can't read frame file \"%s\"\n",str);
		continue;		/* try next one */
            }
            if (first) {
                strcpy(firstfile,str);
                first = 0;
            }
            movie[frame] = pr_load(fp,cmap);
            fclose(fp);
            fprintf(stderr,"rasterfile %d: %s\n",frame,str);
            frame++;
        }
        frame++; /* overshoot like above */
    }

    fprintf(stderr,"\n");
    nframe = frame - 1;
    lastframe = frame - 2;
    frame = 0;
    if (nframe == 0)
        fprintf(stderr,"No rasterfiles found. Exiting\n");
    return(nframe);
}

Movie_Fini()
{
    EraseScreen();
    (void)fclose(mouse);
}


Bogus_Movie()
{
    int npx, npy, i, n, frameop;
#if defined(SUNOS4)
    fd_set rd, wd, ed;
#else
    int rd;
#endif
    char c, cmd[81];
    static struct timeval delay = {0, 0};
    int istep;

    base_frame = window_create(NULL, FRAME,
        FRAME_LABEL,      "movie_sv",
        FRAME_SHOW_LABEL, FALSE,
        WIN_ERROR_MSG,    "Can't create window",
        0);

    window_set (base_frame,
        WIN_SHOW,    TRUE,
        WIN_WIDTH,   xsize+10,
        WIN_HEIGHT,  ysize+10,
        0);

    canvas = window_create(base_frame, CANVAS,
        0);
    pw = canvas_pixwin(canvas);
    npx = (int) window_get (canvas, CANVAS_WIDTH);
    npy = (int) window_get (canvas, CANVAS_HEIGHT);

    (void) notify_do_dispatch();        /* implicit dispatch */

    window_set (base_frame, WIN_SHOW, TRUE, 0);     /* show it */
    fflush(stdin);
    printf("MOVIE_SV Demonstration mode, limited command set:\n");
    printf("Type help or ? to get help menu\n");
    fflush(stdout);

    for (i=0; i<=lastframe; i++)         /* init frame labels */
        sprintf(ftext[i],"%d",i);

    i = 0;  /* frame pointer */
    istep = 1;  /* start in forward loop */
    frameop = PIX_SRC;     /* start in normal video (not reverse) */
    for (;;) {           /* LOOP FOREVER */
        pw_write(pw,0,0,npx,npy,frameop,movie[i],0,0);    /* load frame */
        pw_text(pw,5,npy-5,frameop,0,ftext[i]);  /* and id */
#if defined(SUNOS4)
        FD_ZERO(&wd);
        FD_ZERO(&ed);
        FD_ZERO(&rd);
        FD_SET(0,&rd);      /* only check for input */
        if (select(1,&rd, &wd, &ed, &delay)) {
            if (FD_ISSET(0,&rd)) {      /* check if input */
#else        
        rd = 0x01;  /* only check for input */
        if (select(1,&rd,0,0,&delay)) {
            if (rd & 0x01) {        /* check if input */
#endif
                n=read(0,cmd,80);     /* and read it */
                if (n>1) {      /* discard just newlines */
                    cmd[n-1] = '\0';    /* get rid of newline */
                    if (cmd[0]=='q')
                        break;      /* quit */
                    else if (cmd[0]=='?' || cmd[0]=='h') {
                        Bogus_Help();
                        i -= istep;
                    } else if (cmd[0]=='+' || cmd[0]=='=' ||
                               cmd[0]=='.' || cmd[0]=='>' ||
                               strncmp(cmd,"ssf",3)==0) {
                        istep = 1;
                        delay.tv_sec = 10000;
                    } else if (cmd[0]=='-' || cmd[0]=='<' || cmd[0]==',' ||
                                strncmp(cmd,"ssb",3)==0) {
                        istep = -1;
                        delay.tv_sec = 10000;
                    } else if (strncmp(cmd,"ss",2)==0) {        /* SET SPEED */
                        if (strlen(cmd) > 2)
                            n = atoi(&cmd[2]);  /* take argument of ss */
                        else
                            n = 100;            /* take mean */
                        if (n <= -100) n = -99;
                        if (n > 100)  n = 100;
                        framerate = n;
                        delay.tv_sec  = delaytab_s [framerate+100];
                        delay.tv_usec = delaytab_us[framerate+100];
                        printf ("new speed = %d\n",framerate);
                    } else if (strncmp(cmd,"ff",2)==0) {     /* FAST FORWARD */
                        istep = 1;
                        delay.tv_sec = 0;
                        delay.tv_usec = 0;
                    } else if (strncmp(cmd,"fb",2)==0) {    /* FAST BACKWARDS */
                        istep = -1;
                        delay.tv_sec = 0;
                        delay.tv_usec = 0;                        
                    } else if (cmd[0]=='f') {
                        IncreaseRate();
                        delay.tv_sec  = delaytab_s [framerate+100];
                        delay.tv_usec = delaytab_us[framerate+100];
                        printf ("new speed = %d\n",framerate);
                    } else if (cmd[0]=='s') {
                        DecreaseRate();
                        delay.tv_sec  = delaytab_s [framerate+100];
                        delay.tv_usec = delaytab_us[framerate+100];
                        printf ("new speed = %d\n",framerate);
		    } else if (cmd[0]=='r') {		/* REVERSE  VIDEO */
		        if (frameop == PIX_SRC) 
                            frameop = PIX_NOT(PIX_SRC);
                        else 
                            frameop = PIX_SRC;
                    }
                }
            }
            fflush(stdin);
        }
        i += istep;
        if (i > lastframe)
            i=0;                    /* reset */
        if (i < 0)
            i=lastframe;            /* reset */
    }
}

Bogus_Help()
{
    printf ("q,^C      Quit program\n");
    printf ("?,h       This help\n");
    printf ("ssXXX     Set speed at XXX [-100:100]\n");
    printf ("fXXX      Faster by XXX\n");
    printf ("sXXX      Slower by XXX\n");
    printf ("bf        Back and Forth loop\n");
    printf ("ssf + = > single step forwards\n");
    printf ("ssb -   < single step bakwards\n");
    printf ("ff        fast forward\n");
    printf ("fb        fast backward\n");
    printf ("r         reverse video\n");
    printf ("\nIn Demonstration mode some of these commands are not implemented\n");

    fflush(stdout);
}
RefreshMenu()
{
    struct pr_prpos wh;
    int r;
    char pn[64];

    wh.pr = screen;
    wh.pos.x = MenuX1 + 3;
    EraseMenu();
/*
    pr_vector(screen,900,0,900,900,PIX_CLR,0);
    pr_vector(screen,901,0,901,900,PIX_CLR,0);
    pr_vector(screen,902,0,902,900,PIX_CLR,0);
*/
    wh.pos.y = MenuY1 + 14;
    for (r = 0; r < MenuRows; r++) {
        pf_text(wh,PIX_SRC,pfont,menu[page][r]);
        wh.pos.y += MenuHeight;
    }
    pf_text(wh,PIX_SRC,pfont," + ");
    wh.pos.x += PageWidth;
    pf_text(wh,PIX_SRC,pfont," - ");
    wh.pos.x += PageWidth;
    sprintf(pn,"%2d %s",page,menu[page][MenuRows]);
    pf_text(wh,PIX_NOT(PIX_SRC),pfont,pn);
}

EraseMenu()
{
    pr_rop(screen,MenuX1,MenuY1,MenuWidth,MenuY2,PIX_SET,NULLPR,0,0);
}

ReadMouse(x,y,b)
int *x,*y,*b;
{
    char mbuf[4];
    int ms, dx, dy, ol, oc, or, gotany, bits;

    *x = 0;
    *y = 0;
    *b = 0;

    bits = (1 << MouseFD);
    gotany = select(16,&bits,(int *)0,(int *)0,&delay);

    if (bits & (1 << MouseFD)) {
        read(MouseFD,mbuf,3);

        ms = mbuf[0];
        dx = mbuf[1];
        dy = mbuf[2];
        *x = dx;
        *y = dy;
        ol = bl;
        oc = bc;
        or = br;
        bl = 4 & ~ms;
        bc = 2 & ~ms;
        br = 1 & ~ms;
        *b = 0;
        if (bl && !ol) *b = 1;
        if (bc && !oc) *b = 2;
        if (br && !or) *b = 3;
    }
}

MenuHandle(y,b)
int y,b;
{
    int r,oldp;

    r = (y - MenuY1) / MenuHeight;

    if (r == MenuRows) {
        oldp = page;
        if (px > PageP2 && page > 0) page--;
        else if (px < PageP2 && page < (MenuPages - 1)) page++;
        if (page != oldp) {
            DrawMenuCursor(px,my);
            RefreshMenu();
            DrawMenuCursor(px,my);
        }
    }
    else if (b == L) ExecMenu(page,r);
}

ExecMenu(page,row)
int page,row;
{
    (*MenuFunction[page][row])();
}

EraseScreen()
{
    if (screen) pr_rop(screen,0,0,1152,900,PIX_SET,NULLPR,0,0);
}

RedrawEverything()
{
    EraseScreen();
    RefreshMenu();
    RefreshMouseLine();
    DrawMenuCursor(px,my);
}

MenuItem(str)
char *str;
{
    int mr, mp;
    char ms[64];

    sscanf(str+5,"%d %d %[^\n]",&mp,&mr,ms);
    SetMenuItem(mp,mr,ms);
}

SetMenuItem(p,r,str,ptr)
int p,r;
char *str;
int (*ptr)();
{
    char tstr[128];
    int i, l, mc;

    mc = (MenuWidth / 8) - 1;

    strcpy(tstr," ");
    strcat(tstr,str);
    l = strlen(tstr);
    if (l < mc) for (i = l; i < mc; i++) strcat(tstr," ");
    if (l > mc) tstr[mc-1] = ' ';
    strncpy(menu[p][r],tstr,mc);
    MenuFunction[p][r] = ptr;

    if (p == page) {
        DrawMenuCursor(px,my);
        RefreshMenu();
        DrawMenuCursor(px,my);
    }
}

InitMenu()
{
    int p, r, m;
    int BOGUS(), M_quit();
    int DisplayFirstFrame(), DisplayNextFrame();
    int DisplayPrevFrame(), DisplayLastFrame(), SingleFrame();
    int ForwardLoop(), BackwardLoop(), VariableLoop();
    int BackAndForth(), PositionFrame(), ReverseVideo();
    int ForwardPlay(), BackwardPlay(), ScrollFrames();

    for (p = 0; p < MenuPages; p++)
        for (r = 0; r <= MenuRows; r++) {
            strcpy(menu[p][r],"");
            MenuFunction[p][r] = BOGUS;
        }

    m = 0;

    SetMenuItem(0,m++,"First Frame",DisplayFirstFrame);
    SetMenuItem(0,m++,"Forward Frame",DisplayNextFrame);
    SetMenuItem(0,m++,"Backward Frame",DisplayPrevFrame);
    SetMenuItem(0,m++,"Last Frame",DisplayLastFrame);
    SetMenuItem(0,m++,"Single Frame",SingleFrame);
    SetMenuItem(0,m++,"Forward Play Once",ForwardPlay);
    SetMenuItem(0,m++,"Backward Play Once",BackwardPlay);
    SetMenuItem(0,m++,"Forward Loop",ForwardLoop);
    SetMenuItem(0,m++,"Backward Loop",BackwardLoop);
    SetMenuItem(0,m++,"Back and Forth Loop",BackAndForth);
    SetMenuItem(0,m++,"Variable Speed Loop",VariableLoop);
    SetMenuItem(0,m++,"Scroll With Mouse",ScrollFrames);
    SetMenuItem(0,m++,"Position Frame",PositionFrame);
    SetMenuItem(0,m++,"Reverse Frame Video",ReverseVideo);
    SetMenuItem(0,m++,"Quit",M_quit);
}

M_quit()
{
    go = 0;
}

DisplayFirstFrame()
{
    frame = 0;
    DisplayFrame();
}

DisplayNextFrame()
{
    if (frame < lastframe) {
        frame++;
        DisplayFrame();
    }
}

DisplayPrevFrame()
{
    if (frame > 0) {
        frame--;
        DisplayFrame();
    }
}

DisplayLastFrame()
{
    frame = lastframe;
    DisplayFrame();
}

SingleFrame()
{
    int dx, dy, but, Cut(), DisplayPrevFrame(), DisplayNextFrame();

    BindButton(L,"Backward",DisplayPrevFrame);
    BindButton(C,"Stop",Cut);
    BindButton(R,"Forward",DisplayNextFrame);

    rolling = 1;
    while (rolling) {
        ReadMouse(&dx,&dy,&but);
        (*Button[but])();
    }

    BindButton(L,"Select",BOGUS);
    BindButton(C," ",BOGUS);
    BindButton(R," ",BOGUS);
}

ForwardPlay()
{
    if (framerate < 0) framerate *= -1;

    for (frame = 0; frame < nframe; frame++) {
        if (framerate < 0) framerate = 0;
        DisplayFrame();
        DelayFrame();
    }
    frame = lastframe;
}

BackwardPlay()
{
    if (framerate > 0) framerate *= -1;

    for (frame = lastframe; frame >= 0; frame--) {
        if (framerate > 0) framerate = 0;
        DisplayFrame();
        DelayFrame();
    }
    frame = 0;
}

ForwardLoop()
{
    int IncreaseRate(), DecreaseRate(), Cut(), dx, dy, but;

    BindButton(L,"Slower",DecreaseRate);
    BindButton(C,"Stop",Cut);
    BindButton(R,"Faster",IncreaseRate);

    if (framerate < 0) framerate *= -1;

    rolling = 1;
    frame = 0;

    while (rolling) {
        DisplayFrame();
        ReadMouse(&dx,&dy,&but);
        (*Button[but])();
        DelayFrame();
        if (framerate < 0) framerate = 0;
        if (framerate != 0) frame = (frame + 1) % nframe;
    }
    frame = lastframe;
    BindButton(L,"Select",BOGUS);
    BindButton(C," ",BOGUS);
    BindButton(R," ",BOGUS);
}

BackwardLoop()
{
    int IncreaseRate(), DecreaseRate(), Cut(), dx, dy, but;

    BindButton(L,"Slower",DecreaseRate);
    BindButton(C,"Stop",Cut);
    BindButton(R,"Faster",IncreaseRate);

    if (framerate > 0) framerate *= -1;

    rolling = 1;
    frame = 0;

    while (rolling) {
        DisplayFrame();
        ReadMouse(&dx,&dy,&but);
        (*Button[but])();
        DelayFrame();
        if (framerate > 0) framerate = 0;
        if (framerate != 0) frame = (frame + lastframe) % nframe;
    }
    frame = lastframe;
    BindButton(L,"Select",BOGUS);
    BindButton(C," ",BOGUS);
    BindButton(R," ",BOGUS);
}

VariableLoop()
{
    int IncreaseRate(), DecreaseRate(), Cut(), dx, dy, but;

    BindButton(L,"Slower",DecreaseRate);
    BindButton(C,"Stop",Cut);
    BindButton(R,"Faster",IncreaseRate);

    rolling = 1;
    frame = 0;

    while (rolling) {
        DisplayFrame();
        ReadMouse(&dx,&dy,&but);
        framerate += dy;
        if (framerate >  100) framerate =  100;
        if (framerate < -100) framerate = -100;
        (*Button[but])();
        DelayFrame();
        if (framerate > 0) frame = (frame + 1) % nframe;
        else if (framerate < 0) frame = (frame + lastframe) % nframe;
    }


    BindButton(L,"Select",BOGUS);
    BindButton(C," ",BOGUS);
    BindButton(R," ",BOGUS);
}

BackAndForth()
{
    int IncreaseRate(), DecreaseRate(), Cut(), dx, dy, but, dir;

    BindButton(L,"Slower",DecreaseRate);
    BindButton(C,"Stop",Cut);
    BindButton(R,"Faster",IncreaseRate);

    if (framerate < 0) framerate *= -1;

    rolling = 1;
    frame = 0;
    dir = 0;

    while (rolling) {
        DisplayFrame();
        ReadMouse(&dx,&dy,&but);
        (*Button[but])();
        DelayFrame();
        if (framerate < 0) framerate = 0;
        if (framerate != 0) {
            if (dir) {
                if (frame == lastframe) {
                    dir = !dir;
                    frame--;
                }
                else frame++;
            }
            else {
                if (frame == 0) {
                    dir = !dir;
                    frame++;
                }
                else frame--;
            }
        }
    }

    BindButton(L,"Select",BOGUS);
    BindButton(C," ",BOGUS);
    BindButton(R," ",BOGUS);
}

ScrollFrames()
{
    int dx, dy, but, x, Cut();

    BindButton(L,"Stop",Cut);

    x = frame * 8;
    rolling = 1;

    while (rolling) {
        DisplayFrame();
        ReadMouse(&dx,&dy,&but);
        (*Button[but])();

        x += dx;
        frame = x / 8;
        if (frame > lastframe) {
            frame = lastframe;
            x = frame * 8;
        }
        if (frame < 0) {
            frame = 0;
            x = 0;
        }
    }
    BindButton(L,"Select",BOGUS);
}

IncreaseRate()
{
    if (framerate < 100) framerate += 1;
}

DecreaseRate()
{
    if (framerate > -100) framerate -= 1;
}

Cut()
{
    rolling = 0;
}

DisplayFrame()
{
    Annotation();

    if (frame == screenframe) return(0);

    screenframe = frame;
    pr_rop(screen,framex,framey,
        movie[frame]->pr_size.x,movie[frame]->pr_size.y,
        frameop,movie[frame],0,0);
}

Annotation()
{
    char str[256];

    sprintf(str,"Frame %2d of %d   %4d",frame,lastframe,framerate);
    pf_text(xywh,PIX_NOT(PIX_SRC),pfont,str);
}

DelayFrame()
{
    int f, bits, gotsome;
    char mbuf[3];

    bits = 1 << MouseFD;

    f = framerate + 100;
    delay.tv_sec  = delaytab_s[f];
    delay.tv_usec = delaytab_us[f];
    gotsome = select(32,&bits,(int *)0,(int *)0,&delay);
    if (gotsome < 0) {
        perror("select");
        fprintf(stderr,"returned value was %d\n",gotsome);
        DieHorribly();
    }
    if (bits & (1 << MouseFD)) {
        read(MouseFD,mbuf,3);
        framerate += mbuf[2];
        if (framerate > 100) framerate = 100;
        if (framerate < -100) framerate = -100;
        delay.tv_sec = 0;
        delay.tv_usec /= 2;
        select(32,(int *)0,(int *)0,(int *)0,&delay);
    }
}

DieHorribly()
{
    char str[128];

    sprintf(str,"     Goodbye Cruel World!     ");
    pf_text(xywh,PIX_SRC,pfont,str);

    exit(1);
}

PositionFrame()
{
    int dx, dy, but, Cut(), maxx, maxy, oldx, oldy, xs, ys;

    BindButton(L,"Stop",Cut);
    rolling = 1;

    xs = movie[frame]->pr_size.x;
    ys = movie[frame]->pr_size.y;

    maxx = 899 - xs;
    maxy = 899 - ys;

    DrawBox(framex,framey,xs,ys);

    while (rolling) {
        ReadMouse(&dx,&dy,&but);
        (*Button[but])();
        oldx = framex;
        oldy = framey;
        framex += dx;
        framey -= dy;
        if (framex < 0) framex = 0;
        if (framey < 0) framey = 0;
        if (framex > maxx) framex = maxx;
        if (framey > maxy) framey = maxy;
        if (oldx != framex || oldy != framey) {
            DrawBox(oldx,oldy,xs,ys);
            DrawBox(framex,framey,xs,ys);
        }
    }
    EraseFrame();
    DisplayFrame();
    BindButton(L,"Select",BOGUS);
}

ReverseVideo()
{
    if (frameop == PIX_SRC) frameop = PIX_NOT(PIX_SRC);
    else frameop = PIX_SRC;
    EraseFrame();
    DisplayFrame();
}

EraseFrame()
{
    screenframe = -1;
    pr_rop(screen,0,0,900,900,PIX_SET,NULLPR,0,0);
}

DrawBox(x,y,w,h)
int x,y,w,h;
{
    int x2, y2;

    x2 = x + w - 1;
    y2 = y + h - 1;

    pr_vector(screen,x,y,x2-1,y,PIX_NOT(PIX_DST),0);
    pr_vector(screen,x2,y,x2,y2-1,PIX_NOT(PIX_DST),0);
    pr_vector(screen,x2,y2,x+1,y2,PIX_NOT(PIX_DST),0);
    pr_vector(screen,x,y2,x,y+1,PIX_NOT(PIX_DST),0);
}

Touch(f,i)
int f,i;
{
    int pix;

    if (i > 0) {
        if (f < lastframe) pix = pr_get(movie[f+i],0,0);
        else pix = pr_get(movie[0],0,0);
    }
    else {
        if (f > 0) pix = pr_get(movie[f+i],0,0);
        else pix = pr_get(movie[lastframe],0,0);
    }
}
