      blockdata blkdat
c block data segment written by Richard James for 3-D Poisson solver
c
c     initialises labelled common areas
c
      include 'rjinclude.h'
c
      data bdopt,formul,option / 1,2,1 /
      data nspio / 16 /
      data bufer / 40 * 0 /
      data itrbuf / 7 * 0 /
      data iwin, iwout, ibltab / 3, 4, 22 /
      data limfil / 7 /
      data infile, ofile / 14 * -1 /
      data insel, outsel / 14 * 0 /
      data lenn / 0, 0, 0, 0, 0, 0, 0 /
      data iend / 0 /
      data nverif / 12 /
      data iverif / 0, 0, 0, -1, -1, -1,
     1  0, 27, 13, -1, 27, 13,
     2  5,  0, 22,  5, -1, 22,
     3 12, 10,  0, 12, 10, -1,
     4  9, 17, 22,  2,  5, 16,
     5  5,  7,  9,  9, 22, 11,
     6  114 * 0 /
      data smooth / 0.0 /
      data pagsiz / 'small' /
c
c new line added 10 Mar 92
      data pl3 / 1 /
      end
