C
C movie.f for PGPLOT
C
      SUBROUTINE MOVIE0(XMAX)
      REAL*8 XMAX
*
      REAL*4 XMX

      XMX = XMAX
      CALL PGBEGIN(0,'/xs',1,1)
      CALL PGWINDOW(-XMX,XMX,-XMX,XMX)
      CALL PGBOX('bc',0.0,0,'bc',0.0,0)

      RETURN
      END
C
      SUBROUTINE MOVIE1
      CALL PGEND
      RETURN
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
*      CALL plot_all(3, xx(1,1), colour(1),ix,iy,isize)
*
*      write(str,'(F12.4)') TIME
*       Write str to the lower right corner.
*      CALL write_text(str)

*       Flush X buffer.
*      CALL xflush
*

          CALL PGSCI(0)
          CALL PGPOINT(1,SNGL(XX(1,I1)),SNGL(XX(2,I1)),18)
          CALL PGPOINT(1,SNGL(XX(1,I2)),SNGL(XX(2,I2)),18)
          CALL PGPOINT(1,SNGL(XX(1,I3)),SNGL(XX(2,I3)),18)
          CALL PGSCI(2)
          CALL PGSCI(6)
          CALL PGPOINT(1,SNGL(XX(1,I1)),SNGL(XX(2,I1)),1)
          CALL PGSCI(4)
          CALL PGPOINT(1,SNGL(XX(1,I2)),SNGL(XX(2,I2)),1)
          CALL PGSCI(10)
          CALL PGPOINT(1,SNGL(XX(1,I3)),SNGL(XX(2,I3)),1)
          CALL TRANSF(3)
C
C   supposed to copy X into XX, but we don't have X available here yet;
C   see triple.f where this was supposed to happen
C
C
          CALL PGSCI(2)
          CALL PGPOINT(1,SNGL(XX(1,I1)),SNGL(XX(2,I1)),18)
          CALL PGSCI(4)
          CALL PGPOINT(1,SNGL(XX(1,I2)),SNGL(XX(2,I2)),18)
          CALL PGSCI(10)
          CALL PGPOINT(1,SNGL(XX(1,I3)),SNGL(XX(2,I3)),18)
	

      END
