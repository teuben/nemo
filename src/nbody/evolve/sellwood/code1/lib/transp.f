      subroutine transp(x, nsets, ncont)
c routine written by Richard James for 3-D Poisson solver
      integer nsets, ncont
      real x( nsets * ncont )
c
c     purpose
c
c        to transpose a 2-dimensional matrix held in a 1-dimensional
c        array.
c
c     parameters
c
c     x      -  a real array holding the matrix for transposition
c
c     nsets  =  the number of rows in the matrix, each row occupying
c               a contiguous set of locations.  this terminology
c               deviates from normal fortran terminology.
c
c     ncont  =  the length of a contiguous row.
c
c     notes:
c
c     1)   the routine is designed to cope efficiently with extremely
c     elongated matrices.  the elongation may be in the store direction
c     (long row case) or the column direction (long column case).  the
c     path taken in the code differs in these two cases.
c
c     2)   the normal application is to 3-dimensional sets of values,
c     with element (i, j, k) (numbering from zero) having an offset
c     (i*n2*n3 + j*n3 + k) in the long row case, or (j*n1*n2 + k*n1 + i)
c     in the long column case.  entry with nsets = n1, ncont = n2*n3
c     changes the data from the first to the second ordering.  entry
c     with nsets = n2*n3, ncont = n1 reverses this operation.
c
      include 'rjinclude.h'
c
c local variables
      integer i, is, iss, ivec, nelem
c
c      real x(nsets*ncont)
c
c     calculate number of elements and check stack size.
c
      nelem = nsets*ncont
      i = istack + nelem
      if(i.gt.maxstk) then
        maxstk = i
        mxstid = 'transp'
      end if
      if(i.gt.lstack) then
        write( s2, 200 )istack, nelem, lstack
 200    format( ' stack overflow in transp(xa)' /
     1          ' stack pointer =', i10, ' number of values =',
     2          i10, ' limit =', i10 )
        call crash( 'transp', 'xa 1')
        stop
      end if
c
c     decide if given matrix is long row or long column.
c
      if(ncont.lt.nsets) go to 8
c
c     copy matrix to work area.
c
      do 1 i = 1, nelem
 1    wstack(istack + i) = x(i)
c
c     transpose long row matrix.
c
      is = 0
      iss = istack - ncont
      do 4 ivec = 1, nsets
      is = is + 1
      iss = iss + ncont
      do 4 i = 1, ncont
 4    x(is + i*nsets - nsets) = wstack(iss + i)
c
c     finish
c
      return
c
c     long column matrix.
c     cycle over vectors for collection
c
 8    is = 0
      iss = istack - nsets
      do 9 ivec = 1, ncont
c
c     pick up vectors
c
      is = is + 1
      iss = iss + nsets
      do 9 i = 1, nsets
 9    wstack(iss + i) = x(is + i*ncont - ncont)
c
c     copy back to matrix and finish
c
      do 7 i = 1, nelem
 7    x(i) = wstack(istack + i)
      return
      end
