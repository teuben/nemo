#include <pixrect/pixrect_hs.h>
#include <stdio.h>
#include <sgtty.h>
#include <fcntl.h>
#include <sys/time.h>

#define NULLPR (struct pixrect *)NULL

struct pixrect *screen;
struct pixfont *pf_default(), *pfont;

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
#define MAXFRAMES 512
struct pixrect *movie[MAXFRAMES];
int nframe, lastframe, frame, rolling, screenframe, framerate;
struct timeval delay;
int delaytab_s[201], delaytab_us[201];
int framex, framey, frameop;

char menu[MenuPages][MenuRows+1][128], MouseString[4][11];
int MainX1, MainY1, MainX2, MainY2;

char scriptfile[MAXFRAMES];

main(argc,argv)
int argc;
char *argv[];
{
    int argn;

    scriptfile[0] == '\0';

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
                default:
                    Usage();
                    exit(0);
            }
        }
        else {
            Usage();
            exit(0);
        }
    }

    Init();
    Menu();
    Fini();
}

Usage()
{
    fprintf(stderr,"movie [-p %%d %%d] [-s %%s]\n");
    fprintf(stderr,"    -p   Frame row col position\n");
    fprintf(stderr,"    -s   Script file name\n");
}

Menu()
{
    int dx, dy, but;

    DrawMenuCursor(px,my);

    while (go) {
        ReadMouse(&dx,&dy,&but);
        MoveMenuCursor(dx,dy);
        if (but) MenuHandle(my,but);
    }

    DrawMenuCursor(px,my);
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

Init()
{
    FILE *fopen();
    int one = 1;
    int BOGUS(), Menu();
    int i, s, us;
    double d, dd, ds, dus;

    LoadFrames();
    screen = pr_open("/dev/fb");
    mouse = fopen("/dev/mouse","r");
    if (mouse == NULL) {
        (void)printf("can't open /dev/mouse\n");
        exit(1);
    }
    MouseFD = fileno(mouse);

    go = 1;
    framerate = 90;
    screenframe = -1;
    frame = 0;
    framex = 0;
    framey = 0;
    frameop = PIX_NOT(PIX_SRC);

    delaytab_s [  0] = delaytab_s [200] = 0;
    delaytab_us[  0] = delaytab_us[200] = 0;
    d = 0.75;
    dd = d  / 99.0;

    for (i = 1; i < 100; i++) {
        s = (int)d;
        dus = (d - (double)s) * 100000.0;
        us = (int)dus;
        delaytab_s[i + 100] = delaytab_s[100 - i] = s;
        delaytab_us[i + 100] = delaytab_us[100 - i] = us;
        d -= dd;
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
}

LoadFrames()
{

    int ifile;
    char str[256], *malloc(), *cmap;
    FILE *fp, *fopen(), *script;
    struct pixrect *pr_load();
/*    struct colormap_t *cmap; */

    fprintf(stderr,"Loading frames\n");
/*    cmap = (struct colormap_t *) malloc( sizeof(struct colormap_t) ); /* PJT */

    cmap = malloc(5*4); /* kludge */




    if (scriptfile[0] == '\0') {
        ifile = 1;

        for (frame = 0; frame < MAXFRAMES && ifile; frame++) {
            sprintf(str,"frame.%d",frame);
            fp = fopen(str,"r");
            if (fp == NULL) ifile = 0;
            else {
                movie[frame] = pr_load(fp,cmap);
                fclose(fp);
                fprintf(stderr,"%2d ",frame);
                if (9 == (frame % 10)) fprintf(stderr,"\n");
                fflush(stderr);
            }
        }
    }
    else {
        script = fopen(scriptfile,"r");
        if (script == NULL) {
            fprintf(stderr,"movie: can't read script file\"%s\"\n",scriptfile);
            exit(1);
        }

        frame = 0;
        while (EOF != fscanf(script,"%s",str)) {
            fp = fopen(str,"r");
            if (fp == NULL) {
                fprintf(stderr,"movie: can't read frame file \"%s\"\n",str);
                exit(1);
            }
            movie[frame] = pr_load(fp,cmap);
            fclose(fp);
            fprintf(stderr,"frame %d: %s\n",frame,str);
            frame++;
        }
        frame++; /* overshoot like above */
    }

    fprintf(stderr,"\n");
    nframe = frame - 1;
    lastframe = frame - 2;
    frame = 0;
    if (nframe == 0) Fini();
}

Fini()
{
    EraseScreen();
    (void)fclose(mouse);
    exit(0);
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
    pr_rop(screen,0,0,1152,900,PIX_SET,NULLPR,0,0);
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
