/*
 * SNAP3DV:  convert snapshot to 3dv format for 3D display
 *           (the output format is described in the manual page)
 *
 *      26-sep-91   Created             Peter Teuben
 *	 9-oct-91   full NEMO program with snapshots         PJT
 *	10-oct-91   option to write out AcroSpin(TM) format  PJT
 *	 3-feb-95   1.0c: simple xyzc format		     PJT
 *       7-sep-95   1.0d: added x3d format for point sources PJT
 *	20-sep-00   1.0e: added speck format for partiview   PJT
 */

#include <stdinc.h>
#include <getparam.h>
#include <vectmath.h>		/* otherwise NDIM undefined */
#include <filestruct.h>

#include <snapshot/snapshot.h>	
#include <snapshot/body.h>
#include <snapshot/get_snap.c>


string defv[] = {
    "in=???\n           Input snapshot",
    "out=???\n          Output ascii file for 3DV",
    "xvar=x\n           X-bodytrans variable to view",
    "yvar=y\n           Y-bodytrans variable to view",
    "zvar=z\n           Z-bodytrans variable to view",
    "xrange=\n          X-range of variables [autoscale]",
    "yrange=\n          Y-range of variables [autoscale]",
    "zrange=\n          V-range of variables [autoscale]",
    "time=\n		Time to select (default: first one)",
    "visib=1\n          Visibility bodytrans (0=no >=yes & layer number)",
    "color=3\n          Color of particles (bodytrans)",
    "border=0\n         Color and Color-1 of borders (0=not plotted)",
    "format=%f\n        Format used to store floating point coords X,Y,Z",
    "mode=3dv\n         Output mode {3dv, acd, dxf, wld, xyzc, x3d}",
    "VERSION=1.0c\n     3-feb-92 PJT",
    NULL,
};

string usage="convert snapshot for 3d viewing formats";

local void setrange(real *, string);

void nemo_main(void)
{
    stream instr, fo;
    char fmt[80];
    int i, n, color, nbody, bits, visib, nout, border, major_col, minor_col;
    real x,y,z, xmin,ymin,zmin, xmax,ymax,zmax, tsnap, tselect;
    real xrange[2], yrange[2], zrange[2];
    string xvar, yvar, zvar, realfmt, timestr, mode;
    rproc btrtrans(), xfunc, yfunc, zfunc;
    iproc btitrans(), cfunc, vfunc;
    Body *btab=NULL, *bp, *bq;

    /* Open files */
    instr = stropen(getparam("in"),"r");
    fo = stropen(getparam("out"),"w");

    /* get various bodytrans functions */
    xvar = getparam("xvar");       xfunc = btrtrans(xvar);
    yvar = getparam("yvar");       yfunc = btrtrans(yvar);
    zvar = getparam("zvar");       zfunc = btrtrans(zvar);
    cfunc = btitrans(getparam("color"));
    vfunc = btitrans(getparam("visib"));

    /* various remaining parameters */
    border = getiparam("border");
    realfmt = getparam("format");
    sprintf(fmt,"%s %s %s",realfmt,realfmt,realfmt);
    dprintf(0,"Format used in ascii 3dv output file=%s\n",realfmt);
    timestr = getparam("time");
    if (timestr!=NULL && *timestr!=NULL) tselect=getdparam("time");
    mode = getparam("mode");

    get_history(instr);
    for(n=0;;) {  /* read until first matched time found ; n = #stars found */
	get_history(instr);                         /* paranoia: */
        if (!get_tag_ok(instr, SnapShotTag))        /* if no more snapshots */
            break;                                  /* done with work */
        get_snap(instr, &btab, &nbody, &tsnap, &bits);  /* get snapshot */
        if(timestr!=NULL && *timestr!=NULL) {       /* if time select */
           if (tsnap != tselect) continue;		/* skip this time */
	}
        for (bp=btab, bq=btab, i=0; bp<btab+nbody; bp++, i++) {  /* select */
            x = xfunc(bp,tsnap,i);
            y = yfunc(bp,tsnap,i);
            z = zfunc(bp,tsnap,i);
            color = cfunc(bp,tsnap,i);      /* color of particle */
            visib = vfunc(bp,tsnap,i);      /* also called: the layer */
            if (visib>0) {                  /* visib can be any > 0 */
                if (n==0) {
                    xmin=xmax=x;
                    ymin=ymax=y;
                    zmin=zmax=z;
                }
                xmin=MIN(x,xmin); xmax=MAX(x,xmax);
                ymin=MIN(y,ymin); ymax=MAX(y,ymax);
                zmin=MIN(z,zmin); zmax=MAX(z,zmax);
                Pos(bq)[0] = x;
                Pos(bq)[1] = y;
                Pos(bq)[2] = z;
                Key(bq) = color;                   /* remember Key as color */
                Aux(bq) = visib;    /* remember Aux (a real) as layer/visib */
                if(bp!=bq)
                    bcopy(bp,bq,sizeof(Body));
                n++;
                bq++;
            } /* if(visib) */
        } /* for (bp) */
        dprintf(0,"Copied %d particles for 3dv\n",n);

        if (1) break;                               /* only do one for now */

    } /* for(;;) */
    strclose(instr);
    if (n==0) { 
	strdelete(fo,TRUE);
	error("Last read snapshot at time=%f with nbody=%d no match",
		tsnap,nbody);
    }

    dprintf(0,"Found min and max:\n  %s=%f %f\n  %s=%f %f\n  %s=%f %f\n",
            xvar,xmin,xmax,yvar,ymin,ymax,zvar,zmin,zmax);
    setrange(xrange,getparam("xrange"));
    setrange(yrange,getparam("yrange"));
    setrange(zrange,getparam("zrange"));
    if(xrange[0]!=xrange[1]) { xmin=xrange[0]; xmax=xrange[1]; }
    if(yrange[0]!=yrange[1]) { ymin=yrange[0]; ymax=yrange[1]; }
    if(zrange[0]!=zrange[1]) { zmin=zrange[0]; zmax=zrange[1]; }

    for(bp=btab, bq=btab, nout=0; bp<btab+n; bp++) {
        if(xmin>Pos(bp)[0] || Pos(bp)[0]>xmax ||
           ymin>Pos(bp)[1] || Pos(bp)[1]>ymax ||
           zmin>Pos(bp)[2] || Pos(bp)[2]>zmax)   continue;
        if(bp!=bq)
            bcopy(bp,bq,sizeof(Body));
        bq++;
        nout++;
    }
    if(n!=nout) printf("%d particles outside selected box\n",n-nout);
    if(nout==0) warning("No particles left to plot");
    if(nout>1500) warning("Too many particles for a good 3DV run");
    n=nout;

/** Now the different modes will be handled **/
  if (streq(mode,"3dv") || streq(mode,"3DV")) {             /** 3DV format **/


    if(border>0)				/* print header */
        fprintf(fo,"%d\n",n+8);
    else
        fprintf(fo,"%d\n",n);

    for (bp=btab, i=0; bp<btab+n; bp++, i++) {	/* all XYZ of stars */
        fprintf(fo,fmt,Pos(bp)[0],Pos(bp)[1],Pos(bp)[2]);  
        fprintf(fo,"\n");
    }

    if(border>0){				/* all XYZ of bounding box */
        fprintf(fo,fmt,xmin,ymin,zmin); fprintf(fo,"\n");
        fprintf(fo,fmt,xmin,ymin,zmax); fprintf(fo,"\n");
        fprintf(fo,fmt,xmin,ymax,zmin); fprintf(fo,"\n");
        fprintf(fo,fmt,xmax,ymin,zmin); fprintf(fo,"\n");
        fprintf(fo,fmt,xmax,ymax,zmax); fprintf(fo,"\n");
        fprintf(fo,fmt,xmax,ymax,zmin); fprintf(fo,"\n");
        fprintf(fo,fmt,xmax,ymin,zmax); fprintf(fo,"\n");
        fprintf(fo,fmt,xmin,ymax,zmax); fprintf(fo,"\n");
    }
    if (border>0)				/* header of segments */
        fprintf(fo,"%d\n",2*n+16);
    else
        fprintf(fo,"%d\n",2*n);   

    for (bp=btab, i=0 ; i<n; bp++, i++) {    	/* move and draw all stars */
        fprintf(fo,"%d 0\n%d %d\n",i+1,i+1,Key(bp));
    }
    if (border>0) {				/* draw the bounding box */
        major_col = border;             /* kludge */
        minor_col = border+1;
        fprintf(fo,"%d %d\n",n+2,0);
        fprintf(fo,"%d %d\n",n+1,major_col);
        fprintf(fo,"%d %d\n",n+3,major_col);
        fprintf(fo,"%d %d\n",n+1,0);
        fprintf(fo,"%d %d\n",n+4,major_col);
        fprintf(fo,"%d %d\n",n+6,minor_col);
        fprintf(fo,"%d %d\n",n+5,minor_col);
        fprintf(fo,"%d %d\n",n+7,minor_col);
        fprintf(fo,"%d %d\n",n+4,minor_col);
        fprintf(fo,"%d %d\n",n+7,0);
        fprintf(fo,"%d %d\n",n+2,minor_col);
        fprintf(fo,"%d %d\n",n+8,minor_col);
        fprintf(fo,"%d %d\n",n+3,minor_col);
        fprintf(fo,"%d %d\n",n+6,minor_col);
        fprintf(fo,"%d %d\n",n+8,0);
        fprintf(fo,"%d %d\n",n+5,minor_col);
    }
  } else if (streq(mode,"acd") || streq(mode,"ACD")) {        /** ACROSPIN **/
    fprintf(fo,"Pointlist X Y Z Color Layer\n");
    for (bp=btab; bp<btab+n; bp++) {	            /* all XYZ of stars */
        fprintf(fo,fmt,Pos(bp)[0],Pos(bp)[1],Pos(bp)[2]);
        fprintf(fo," %d %d\n",Key(bp),(int)Aux(bp));     /* color and layer */
    }
    strclose(fo);
  } else if (streq(mode,"xyzc") || streq(mode,"XYZC")) {          /** XYZC **/
    dprintf(0,"X Y Z Color\n");
    for (bp=btab; bp<btab+n; bp++) {	            /* all XYZ of stars */
        fprintf(fo,fmt,Pos(bp)[0],Pos(bp)[1],Pos(bp)[2]);
        fprintf(fo," %d\n",Key(bp));     		/* color */
    }
    strclose(fo);
  } else if (streq(mode,"x3d") || streq(mode,"X3D")) {            /** X3D **/
    dprintf(0,"NP NS\nX Y Z..\nP1 P2..\n");
    fprintf(fo,"%d %d\n",n,n);
    for (bp=btab; bp<btab+n; bp++) {	            /* all XYZ of stars */
        fprintf(fo,fmt,Pos(bp)[0],Pos(bp)[1],Pos(bp)[2]);
        fprintf(fo,"\n");
    }
    for (i=0; i<n; i++) {	                       /* segments to plot */
        fprintf(fo,"%d %d\n",i,i);
    }
    strclose(fo);
  } else if (streq(mode,"wld") || streq(mode,"WLD")) {         /** VIEWWLD **/
    for (bp=btab; bp<btab+n; bp++) {	            /* all XYZ of stars */
        fprintf(fo,fmt,Pos(bp)[0],Pos(bp)[1],Pos(bp)[2]);
        fprintf(fo," %d\n",2);     /* 2 = 'draw point' command for ViewWld */
    }
    strclose(fo);
  } else if (streq(mode,"dcx") || streq(mode,"DCX")) {         /** AutoCAD **/
    warning("DCX format not supported - use ACD and convert with ACROTRAN");
    strdelete(fo,TRUE);
  } else {
    error("Illegal output mode %s",mode);
    strdelete(fo,TRUE);
  }
}


local void setrange(real *rval, string rexp)
{
    char *cptr;

    if (rexp==NULL || *rexp==0) {
      rval[0] = rval[1] = 0.0;
    } else {
      cptr = strchr(rexp, ':');
      if (cptr != NULL) {
        rval[0] = atof(rexp);
	rval[1] = atof(cptr+1);
      } else {
        rval[0] = 0.0;
	rval[1] = atof(rexp);
      }
    }
}

