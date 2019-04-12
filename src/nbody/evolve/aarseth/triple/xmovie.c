/* 
 * xmovie.c 	C-to-Fortran interface from movie.f to Xtc.c
 *   for Nbody4 
 *   using Xtc.c
 *   Makoto Taiji
 *	linux adaptation by Peter Teuben	jan.98 (Amsterdam)
 *	fixed Y sign convention			23-feb-98	pjt
 *   Fixed second-underscore back to single
 */

/* linux uses double trailing underscores since it has an embedded underscore */
/* or use the g77 compile flag  */
/* -fno-second-underscore	*/

/* gfortran doesn't do this anymore, which is the current default */
#ifdef linux
# define init_movie init_movie_
# define plot_all   plot_all_
# define write_text write_text_
# define x_flush    x_flush_
#else
# define init_movie init_movie__
# define plot_all   plot_all__
# define write_text write_text__
# define x_flush    x_flush_
#endif


static double xmax, ymax, scale;

void
init_movie(size, windowsize)
    double *size;
    int *windowsize;
{
    xmax = ymax = *size;
    scale = *windowsize * 0.5/xmax;
    initgraph("Nbody4",*windowsize,*windowsize);
}

void
plot_all(n, pos, col, ix, iy, size)
    int *n, *col, *ix, *iy, *size;
    double pos[][3];
{
    int i,x,y,ix0, iy0;

    ix0 = *ix-1; iy0 = *iy-1;
    cleardevice_pixmap();
    for(i=0;i<*n;i++,col++) {
      int x,y;
      x = (pos[i][ix0]+xmax)*scale;
      y = (ymax-pos[i][iy0])*scale;
      setcolor_pixmap(*col);
      fillellipse_pixmap(x,y,*size,*size); 
/*      putpixel_pixmap(x,y,*col); */
    }
}

void
write_text(str)
    char * str;
{
    setcolor_pixmap(15);
    outtextxy_pixmap(getmaxx()-100, getmaxy()-20, str);
}

void
x_flush()
{
    copy_from_pixmap();
    xflush();
}
