      subroutine dtrfm3(p, d, d1, d2, synth, cosan)
c routine written by Richard James for 3-D Poisson solver
c
c     purpose
c
c        to change representation for a 2-dimensional mesh along the
c        row direction.  the data may be converted from one
c        repesentation to another, or from fourier to direct form.
c
c     parameters
c
c     p   -  (type integer)  =  start address of data in /scm/
c
c     d   -  (type integer)  =  start address for results in /scm/
c
c     d1  -   (type integer)  =  number of columns in data
c
c     d2  -   (type integer)  =  number of rows in data
c
c     synth - (type logical)  -  indicates that the data must be
c                                synthesised before analysis
c
c     cosan - (type logical)  -  requiests cosine analysis as the
c                                final stage.
c
      include 'rjinclude.h'
c
c local variables
      integer i, id, is
c
      integer p, d, d1, d2
      logical synth, cosan
c
c     check meshes do not overlap unless they are the same
c
      if((p.eq.d).or.(iabs(p - d).ge.(d1*d2))) go to 1
      write(s2, 200) p, d, d1, d2
200   format('0source and destination meshes overlap for subroutine lw'
     1/'0meshes begin at', 2i10, 5x, 'and have dimensions', 2i10)
      stop
 1    if(p.eq.d) go to 2
      is = p + d1*d2 - 1
      id = d
      do 3 i = p, is
      w(id) = w(i)
3     id = id + 1
2     if(synth) go to 4
      if(cosan) call fftcha(d2, d1, d1, d, .true.)
      if(.not.cosan) call fftsna(d2, d1, d1, d, .true.)
      return
4     if(cosan) go to 5
      call fftchs(d2, d1, d1, d, .true.)
      call fftsna(d2, d1, d1, d, .true.)
      return
5     call fftsns(d2, d1, d1, d, .true.)
      call fftcsa(d2, d1, d1, d, .true.)
      return
      end
