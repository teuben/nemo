#include <math.h>
#include "vtc.h"
#define NJMAX (120000)
#define DIM (3)
#define PI M_PI

static int g_cpallet[8][3]={{0,0,0},{0,0,15},{0,15,0},{0,15,15},{15,0,0},
			    {15,0,15},{15,15,0},{15,15,15}};
			  
/**********************************************
	ksort.c 
	  sort program presented by E.Kokubo
***********************************************/

static void
sort_particles(first, last, value, index)
int first, last;
double value[];
int index[];
{
  double ref_value;
  int ref_index;
  int i, j;
  int tmp;

  ref_index = (first + last)/2;
  ref_value = value[index[ref_index]];

  i = first;
  j = last;

  for( ; ; ){
    while(value[index[i]] < ref_value) i++;
    while(value[index[j]] > ref_value) j--;
    if(i >= j) break;
    tmp = index[i];
    index[i] = index[j];
    index[j] = tmp;
    i++;
    j--;
  }
  if (first < i-1) sort_particles(first, i-1, value, index);
  if (j+1 < last) sort_particles(j+1, last, value, index);
}

void initial_animation()
{
    int gd,gm;

/*	set_window_size(650,850);
	set_window_size_pixmap(650,650);*/
    /* !!! */
#if 1
    set_window_size(500,500);
    set_window_size_pixmap(500,500);
#else
    set_window_size(700,700);
    set_window_size_pixmap(700,700);
#endif
    initgraph(&gd, &gm, "");
}

void plot_particle(x,n,c,time)
double x[NJMAX][DIM];
int n;
int c[NJMAX];
double time;
{
    int xmax,ymax;
    int cx,cy;
    char ts[100];
    int izindex[NJMAX];
    double valuez[NJMAX];
    int i;
    int ii,irad=10;
    irad = 10*pow(n/100.0,-1.0/3.0);
    xmax = getmaxx();
    ymax = getmaxy();

    for(i=0;i<n;i++)izindex[i] = i;
    for(i=0;i<n;i++)valuez[i]=x[i][2];
          
    sort_particles(0,n-1,valuez,izindex);

    cleardevice();   
    for(i=0;i<n;i++){

	ii = izindex[n-1-i];
	cx = xmax*x[ii][0]; 
	cy = ymax*x[ii][1];
/*          setcolor(0);
            circle(cx,cy,irad);*/
	if(c[0]==15)setcolor(15);
	if(c[ii]==0)setcolor(14);
	if(c[ii]==1)setcolor(12);
/*          if(x[ii][2]<0.2)fillellipse(cx,cy,5,5);*/
	fillellipse(cx,cy,irad,irad);

    }
    gcvt(time,7,ts);
    setcolor(15);
    circle(xmax/2,ymax/2,xmax/2);
    outtextxy(xmax-50,ymax-50,ts);
    xflush();
}   

int set_pallet(c_org,c_diff)
int c_org,c_diff;
{
    int retnum;
	
    if(c_org!=0){
	retnum=((15-c_diff*1.3) * (c_org)/15);
/*	    printf("c_org %d c_diff %d ret %d\n",c_org,c_diff,retnum);*/
	return retnum;
    }
    return 0;
}


void set_pallet_down(gcnum)
int gcnum;
{
    int i;

    for(i=1;i<=7;i++){
	setrgbpalette(i,set_pallet(g_cpallet[gcnum][0],7-i),
		      set_pallet(g_cpallet[gcnum][1],7-i),
		      set_pallet(g_cpallet[gcnum][2],7-i));
/*	    printf("%d cnum %d %d %d\n",i,set_pallet(g_cpallet[gcnum][0],7-i),
	    set_pallet(g_cpallet[gcnum][1],7-i),set_pallet(g_cpallet[gcnum][2],7-i));*/
    }
}


void set_pallet_up(gcnum)
int gcnum;    
{
    int i;

    for(i=8;i<=14;i++){
	setrgbpalette(i,set_pallet(g_cpallet[gcnum][0],14-i),
		      set_pallet(g_cpallet[gcnum][1],14-i),
		      set_pallet(g_cpallet[gcnum][2],14-i));
    }
}

int ret_cnum(flag,z)
int flag;
double z;
{
    if(flag==0){
	if(z>0.9)      return 1;
	else if(z>0.8) return 2;
	else if(z>0.7) return 3;
	else if(z>0.6) return 4;
	else if(z>0.5) return 5;
	else if(z>0.4) return 6;
	else           return 7;
    }
    else{
	if(z>0.9)      return 8;
	else if(z>0.8) return 9;
	else if(z>0.7) return 10;
	else if(z>0.6) return 11;
	else if(z>0.5) return 12;
	else if(z>0.4) return 13;
	else           return 14;
    }	
}


void plot_particle3D(x,n,c,time)
double x[NJMAX][DIM];
int n;
int c[NJMAX];
double time;
{
    static int ini=0;
    static int xmax1,ymax,xmax2;
    static int xmax1h,ymaxh,xmax2h,xpos_l,xpos_r;
    int cx,cy;
    static double theta,sintheta,costheta;
    double fx,fy,fz;
    char ts[100];
    int izindex[NJMAX];
    double valuez[NJMAX];
    int i;
    static int tposx,tposy;
    int ii;
    static int irad=10;

    if(ini==0){
	ini=1;
	cleardevice();
	set_pallet_down(6);
	set_pallet_up(4);
	irad = 4*pow(n/100.0,-1.0/3.0);
	if(irad==0)irad=1;
	xmax2 = getmaxx_pixmap();
	xmax1 = xmax2/2;
	ymax = getmaxy_pixmap();
	xmax1h = xmax1/2;
	ymaxh = ymax/2;
	  
/*	  theta = 3.0;*/
	theta = 0.0;
	sintheta = sin(theta*PI/180.0);
	costheta = cos(theta*PI/180.0);
	outtextxy(0,getmaxy_pixmap()+40 ,"  ----------------- Plasma Ball --------------------");
	outtextxy(0,getmaxy_pixmap()+80 ,"  by Toshiyuki Fukushige, University of Tokyo, Japan");
	outtextxy(0,getmaxy_pixmap()+120 ,"  to change stereo view type, push right mouse button");

	/* for parallel view of stereogram */
	xpos_r=xmax1;
	xpos_l=0;
	tposx=getmaxx()-200;
	tposy=getmaxy()-40;
	setcolor(15);
	strcpy(ts,"parallel view");
	outtextxy(tposx,tposy,ts);
    }
    cleardevice_pixmap();

    for(i=0;i<n;i++)izindex[i] = i;
    for(i=0;i<n;i++)valuez[i]=x[i][2];
          
    sort_particles(0,n-1,valuez,izindex);
	  
    for(i=0;i<n;i++){
	ii = izindex[n-1-i];
	cx = xmax1*x[ii][0]; 
	cy = ymax*x[ii][1];
	setcolor_pixmap(ret_cnum(c[ii],x[ii][2]));
	    
/*            if(c[ii]==0)setcolor_pixmap(14);
	      if(c[ii]==1)setcolor_pixmap(12);
	      if(x[ii][2]<0.2)fillellipse(cx,cy,5,5);*/

	fx = costheta*(x[ii][0]-0.5)-sintheta*(x[ii][2]-0.5)+0.5; 
	fz = (costheta*(x[ii][2]-0.5)+sintheta*(x[ii][0]-0.5)+0.5)*0.167;
	cx = (xmax1*fx-xmax1h)/(1.0+fz)+xmax1h+xpos_l;
	cy = (ymax*x[ii][1]-ymaxh)/(1.0+fz)+ymaxh;
	fillellipse_pixmap(cx,cy,irad,irad);

	fx = costheta*(x[ii][0]-0.5)+sintheta*(x[ii][2]-0.5)+0.5;
	fz = (costheta*(x[ii][2]-0.5)-sintheta*(x[ii][0]-0.5)+0.5)*0.167;
	cx = (xmax1*fx-xmax1h)/(1.0+fz)+xmax1h+xpos_r;
	cy = (ymax*x[ii][1]-ymaxh)/(1.0+fz)+ymaxh;
	fillellipse_pixmap(cx,cy,irad,irad);
    }
/*          gcvt(time,7,ts);*/
    sprintf(ts,"%6.1f",time);
    setcolor_pixmap(15);
    circle_pixmap(xmax1h,ymaxh,xmax1h);
    circle_pixmap(xmax1h+xmax1,ymaxh,xmax1h);
    outtextxy_pixmap(xmax2-90,ymax-30,ts);

    copy_from_pixmap(0,0);
	 
    if(xgetbutton_now()==3){
	setcolor(0);
	strcpy(ts,"parallel view   ");
	bar(tposx,tposy,tposx+20*textheight(ts),tposy+textwidth(ts));
	setcolor(15);
	if(xpos_r==0){
	    strcpy(ts,"parallel view");
	    outtextxy(tposx,tposy,ts);
	    xpos_r=xmax1;
	    xpos_l=0;
	}
	else{
	    strcpy(ts,"crossing view");
	    outtextxy(tposx,tposy,ts);
	    xpos_l=xmax1;
	    xpos_r=0;
	}
    }
		  
    xflush();
}   

void plot_particle2D(x,n,c,time)
double x[NJMAX][DIM];
int n;
int c[NJMAX];
double time;
{
    static double xtmp[NJMAX][DIM];
    static int ini=0;
    static int xmax1,ymax,xmax2;
    static int xmax1h,ymaxh,xmax2h,xpos_l,xpos_r;
    int cx,cy;
    static double theta,sintheta,costheta;
    double fx,fy,fz;
    char ts[100];
    static int izindex[NJMAX];
    static double valuez[NJMAX];
    int i,d;
    static int tposx,tposy;
    int ii;
    static int irad=10;

    for(i=0;i<n;i++){
	for(d=0;d<DIM;d++) xtmp[i][d] = 1.0*x[i][d];
    }	  

    if(ini==0){
	ini=1;
	initial_animation();
	cleardevice();
	set_pallet_down(6);
	set_pallet_up(4);
	irad = 2.0*pow(n/100.0,-1.0/3.0);
	if(irad==0)irad=1;
	xmax1 = getmaxx_pixmap();
	ymax = getmaxy_pixmap();
	xmax1h = xmax1/2;
	ymaxh = ymax/2;
	  
/*	  theta = 3.0;*/
	theta = 30.0;
	sintheta = sin(theta*PI/180.0);
	costheta = cos(theta*PI/180.0);
    }
    cleardevice_pixmap();

    for(i=0;i<n;i++)izindex[i] = i;
    for(i=0;i<n;i++)valuez[i]=xtmp[i][2];
          
    sort_particles(0,n-1,valuez,izindex);
	  
    for(i=0;i<n;i++){
	ii = izindex[n-1-i];
	cx = xmax1*xtmp[ii][0]; 
	cy = ymax*xtmp[ii][1];
	setcolor_pixmap(ret_cnum(c[ii],xtmp[ii][2]));
	fx = costheta*(xtmp[ii][0]-0.5)-sintheta*(xtmp[ii][2]-0.5)+0.5; 
	fz = (costheta*(xtmp[ii][2]-0.5)+sintheta*(xtmp[ii][0]-0.5)+0.5)*0.167;
	cx = (xmax1*fx-xmax1h)/(1.0+fz)+xmax1h+xpos_l;
	cy = (ymax*xtmp[ii][1]-ymaxh)/(1.0+fz)+ymaxh;
	if (ii < 3) {
	    setcolor_pixmap(14);
	    fillellipse_pixmap(cx,cy,10,10);
	}
	else {
	    fillellipse_pixmap(cx,cy,irad,irad);
	}
    }
    sprintf(ts,"%6.1f",time);
    setcolor_pixmap(15);
    outtextxy_pixmap(xmax1-120,ymax-30,ts);
    copy_from_pixmap(0,0);
    xflush();
}   

void plot_star(x,n,time,ratio,m,initm)
double (*x)[DIM];
int n;   
double time;
double ratio;
double *m,initm;
{
    static double anix[NJMAX][DIM],maxx,button;
    static int tmpc[NJMAX];
    int i,d;
    int num = 0;

    for(i=0;i<n;i++){
#if 0
	num = n;
	for(d=0;d<DIM;d++) anix[i][d] = x[i][d]* ratio + 0.5;
	tmpc[i] = 0;  

	if(m[i]>initm) tmpc[i]=14;
	if(m[i]>2*initm) tmpc[i]=13;
	if(m[i]>0.1) tmpc[i]=12;
	if(m[i]>0.2) tmpc[i]=11;
	if(m[i]>0.3) tmpc[i]=10;    
#else
	if(m[i]<=initm) {
	    anix[num][0] = x[i][0]* ratio + 0.5;
	    anix[num][1] = x[i][1]* ratio + 0.5;
	    anix[num][2] = x[i][2]* ratio + 0.5;
	    tmpc[num] = 0;  
	    num++;
	}
#endif
	/*	printf("tmpc %d initm %g m %g\n",tmpc[i],initm,m[i]);*/
    }
    plot_particle2D(anix,num,tmpc,time);
}       
