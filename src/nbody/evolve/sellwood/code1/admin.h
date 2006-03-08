c
c include file:  admin.h
c
      integer chkstp, ilist, inext, istep, izone, lpf, mlist, mesh
      integer mzone, nbod, ncoor, ncheck, ndistf, ngx, ngxy, ngy, ngz
      integer ni, nlists, no, noff, nplanes, nres, nwpp, outstp
      integer ipt( 8 ), islist( 2, 257 )
      real claus, c3pmf, ke, lscale, pe, pmass, tmass, ts, xm, ym, zm
      real amom( 3 ), com( 3 ), dh( 3 ), lmom( 3 )
      common / admin / ngx, ngy, ngz, ngxy, dh, mesh,
     +                 xm, ym, zm, nplanes, c3pmf,
     +                 izone, inext, ilist, mzone, mlist,
     +                 ni, no, ndistf, nres, outstp,
     +                 lscale, ts, nbod, pmass, tmass, istep,
     +                 chkstp, ncheck, claus, ke, pe, amom, com, lmom,
     +                 ncoor, lpf, nwpp, noff, nlists, ipt, islist
c
      integer lptcls, iptcls( 1 )
      real ptcls( 1 )
      common / ptcls / lptcls, ptcls
      equivalence ( ptcls( 1 ), iptcls( 1 ) )
c
      real coor( 1 )
      common / results / coor
c
      integer lgrids
      real grids( 1 )
      common / grids / lgrids, grids
c
      integer mbuff
      parameter ( mbuff = 100000 )
      integer iflag( mbuff ), iz( mbuff ), loc( mbuff ), ncl( mbuff )
      logical nskip( mbuff )
      real acc( 3, mbuff ), gpot( mbuff ), newc( 6, mbuff )
      real oldc( 6, mbuff ), wt( 8, mbuff )
      common / buffer / oldc, newc, ncl, wt, acc,
     +                  iz, loc, iflag, gpot, nskip
