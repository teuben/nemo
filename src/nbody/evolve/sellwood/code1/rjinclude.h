c
c include file:  rjinclude.h
c
      integer rjn1, rjn2, rjn3
c     parameter ( rjn1 = 257, rjn2 = 257, rjn3 = 257 )
c     parameter ( rjn1 = 129, rjn2 = 129, rjn3 = 129 )
c     parameter ( rjn1 =  65, rjn2 =  65, rjn3 =  65 )
      parameter ( rjn1 =  33, rjn2 =  33, rjn3 =  33 )
c
c / rjarea /
      integer lnarea
c 8192 is value of  lenbuf
      parameter ( lnarea = 4 * 8192 +
     +                 2 * ( rjn2 * rjn3 + rjn1 * rjn3 + rjn1 * rjn2 ) )
      real area( lnarea )
      common / rjarea / area
c / rjblank /
      integer lstack
      parameter ( lstack = rjn1 * rjn2 * rjn3 + 19 * rjn2 * rjn3
     +                     + 6 * ( rjn1 - 1 ) * ( rjn2 + rjn3 ) )
      integer istack, jstack( lstack ), maxstk
      real wstack( lstack )
      common / rjstack / istack, maxstk, wstack
      equivalence ( jstack, wstack )
c
      character mxstid*6
      common / rjmxstid / mxstid
c / rjchan /
      integer ibltab, iend, iwin, iwout, limb, pll4, pl1, pl2, pl3, pl4,
     +        pp, p1, p2, p3, s7, s8, s9
      logical file, fileb, files, filep, meshpr
      common / rjchan / s7, s8, s9, file, fileb, filep, files, p1, p2,
     +                  p3, limb, pl1, pl2, pl3, pl4, pp, pll4, iend,
     +                  iwin, iwout, ibltab, meshpr
c / rjctrl /
      integer base, bd( 20 ), bdopt, brick2, brick3, by, ctrlc, dimm,
     1        formul, lmc, length, m1, m2, m3, n1, n2, n3, option,
     2        scmbuf, skip, s1, s2, s3, s4, s5, s6, schan( 6 )
      logical twodim
      real wt1, wt2, wt3
      common / rjctrl / option, length, brick2, brick3, skip, base,
     1                  scmbuf, n1, n2, n3, m1, m2, m3, twodim, s1, s2,
     2                  s3, s4, s5, s6, lmc, dimm, wt1, wt2, wt3,
     3                  formul, bdopt, ctrlc, by, bd
      equivalence ( schan( 1 ), s1 )
c / rjfactor /
      integer lnfac
      parameter ( lnfac = 2 * rjn3 + ( rjn3 + 1 ) / 2 )
      real fact( lnfac )
      common / rjfactor / fact
c / rjfnamtb /
      character fnamtb( 100 )*8
      common / rjfnamtb / fnamtb
c / rjgreen /
      integer lngrn
      parameter ( lngrn = rjn2 * ( rjn3 + 1 ) )
      real grein( lngrn )
      common / rjgreen / grein
c / rjgrbuf /
      integer icr3, icr6, inx3, inx6, iptgrn, irec3, irec6, limgrn,
     1        nrec3, nrec6
      logical fil3
      common / rjgrbuf / nrec3, nrec6, irec3, irec6, icr3, inx3, icr6,
     1                   inx6, fil3, iptgrn, limgrn
c / rjgrnbuf /
      integer lngrnb
      parameter ( lngrnb = rjn1 * rjn2 * rjn3 +
     +                 2 * ( rjn2 * rjn3 + rjn1 * rjn3 + rjn1 * rjn2 ) )
      logical grnmem
      real grnbuf( lngrnb )
      common / rjgrnbuf / grnmem, grnbuf
c / rjiocon /
      integer itrbuf( 7 ), limfil, nspio
      integer bufer( 2, 20 ), infile( 7 ), inlim( 7 ), inpt( 7 ),
     1        inrec( 7 ), insel( 7 ), inxt( 7 ), lenn( 7 ), lin( 7 ),
     2        lout( 7 ), ofile( 7 ), outpt( 7 ), outnxt( 7 ),
     3        outsel( 7 )
      common / rjiocon / inpt, outpt, inlim, inrec, insel, outsel, inxt,
     1                   outnxt, infile, ofile, lin, lout, lenn, limfil,
     2                   nspio, itrbuf, bufer
c / rjpagsiz /
      character pagsiz*5
      common / rjpagsiz / pagsiz
c / rjlev3 /
      integer lnlev3
      parameter ( lnlev3 = rjn1 * rjn2 * rjn3 )
      integer levend
      real jfile( 192 ), lev3( lnlev3 )
      common / rjlev3 / jfile, lev3, levend
c / rjmarks /
      integer b( 3 ), b3, edge( 3 ), edpl( 3, 2 ), incrl, limit( 3 ),
     2        line, mark( 3 ), nom( 3 ), norm( 3 ), plane, pl( 3, 4 ),
     3        tran
      logical massto, copy, opt3, gren, permit
      real cf( 6 )
      common / rjmarks / plane, line, tran, b, nom, norm, limit, mark,
     1                   edge, pl, b3, edpl, incrl, massto, copy, opt3,
     2                   gren, permit, cf
c / rjmass /
      real factor, scfact
      common / rjmass / factor, scfact
c / rjptr /
      integer q1, q2, q3, q4, q5, q6, q7, q8, q9, rr4, rr5, rr6, rr7,
     1        rr8, rr9,r1, r2, r3, r4, r5, r6, r7, r8, r9
      real corner( 8 ), w1, w2, w3
      integer ii( 24 )
      common / rjptr / q1, q2, q3, q4, q5, q6, q7, q8, q9, r1, r2, r3,
     +                 r4, r5, r6, r7, r8, r9, rr4, rr5, rr6, rr7, rr8,
     +                 rr9, w1, w2, w3, corner
      equivalence ( ii, q1 )
c / rjscm /
      integer lnscm
      parameter ( lnscm = rjn1 * rjn2 * rjn3 +
     +                 2 * ( rjn2 * rjn3 + rjn1 * rjn3 + rjn1 * rjn2 ) )
      real w( lnscm )
      common / rjscm / w
c / rjsmooth /
      real smooth
      common / rjsmooth / smooth
c / rjvect /
      integer lnvect
      parameter ( lnvect = 2 * rjn2 * rjn3  )
      integer ivect( lnvect )
      logical lvect( lnvect )
      real vect( lnvect )
      common / rjvect / vect
      equivalence ( vect, ivect ), ( vect, lvect )
c / rjverify /
      integer iverif( 3, 50 ), nverif
      logical grdir, verify( 3 )
      real val( 2, 50 )
      common / rjverify / nverif, verify, grdir, iverif, val
c
      integer lnwork
c this is more than the bare minimum required - extra space is
c   used to improve efficiency of boundary calculation
      parameter ( lnwork = 2 * rjn2 * rjn3 + rjn1 * rjn2 )
c bare minimum is     
c      lnwork = 2 * rjn2 * rjn3 + 4 * ( rjn2 - 2 ) + rjn3 - 2 )
      integer iwork( lnwork )
      real wfill( 1 ), work1( rjn3 ), work2( rjn3 ), work3( rjn3 ),
     1     work( lnwork )
      common / rjwork / work1, work2, work3, wfill
c      equivalence ( work( 1 ), work( 1 ) ), ( work1( 1 ), iwork( 1 ) )
      equivalence ( work1( 1 ), iwork( 1 ) )
