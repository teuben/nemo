C
C movie.f for X (Xtc)
C
      SUBROUTINE MOVIE0(XMAX)
      REAL*8  XMAX
*
*
*       Initialization for movie with X11.
*       ----------------------------------
*

      INTEGER iwin
*
*
*       Read movie parameters.
*       xmax = distance from centre to the edge.
*       iwin = real window size (in pixels).
*      READ (5,*)  XMAX, IWIN
	IWIN = 512
*
*       Initialize movie window.
*
      call init_movie(xmax,iwin)
*
      RETURN
*
      END
C
      SUBROUTINE MOVIE1
      STOP
      END
C
      SUBROUTINE MOVIE(I1,I2,I3,XX,TIME)
*
*
*       Real time movie with X11.
*       -------------------------
*
*       Developed by Makoto Taiji (Sept/95).
*       ************************************
*
      REAL*8  XX(3,3),TIME
      INTEGER colour(3)
      CHARACTER*16 str
*
*
*       Colour table definition.
*          0 black
*          1 blue
*          2 green 
*          3 cyan 
*          4 red 
*          5 magenta 
*          6 brown 
*          7 light gray 
*          8 gray 
*          9 light blue 
*         10 light green 
*         11 light cyan 
*         12 light red 
*         13 light magenta  
*         14 yellow 
*         15 white 
*
*       Specify colours according to particle name.
      colour(I1) = 15
      colour(I2) = 12
      colour(I3) = 3
*
*       Set coordinate indices and size of points.
      ix = 1
      iy = 2
      isize = 2
*       colour(I)    = colour of particles.
*       ix           = index for abscissa.
*       iy           = index for ordinate.
*       isize        = size of points.
*
      CALL plot_all(3, xx(1,1), colour(1),ix,iy,isize)
*
      write(str,'(F12.4)') TIME
*       Write str to the lower right corner.
      CALL write_text(str)

*       Flush X buffer.
      CALL x_flush
*
      RETURN
*
      END
