c***********************************************************************
c+      Keywords for NEMO retrieved with the ftoc(1NEMO) program.
c   x0=1\n                  Initial x coordinate
c   y0=0\n                  Initial y coordinate
c   u0=0\n                  Initial u velocity
c   v0=1\n                  Initial v velocity
c   per=3\n                 Period
c   type=1\n                Type of orbit {1,2,3}
c   ome=0.0\n               Pattern speed of potential
c   norbit=1\n              Number of orbits to compute
c   step=0.1\n              Step to increase non-zero position
c   VERSION=0.1\n           23-dec-2023 PJT
c-      
c-----------------------------------------------------------------------
c  HENYEY: Compute periodic orbits using the Henyey method
c          of searching for periodic solutions of a system of 
c          differential equations in a rotating 2D potential.
c
c          The method has been used first by T.S. v. Albada & 
c          R.H. Sanders (1982) MNRAS, 201, 303. 
c          For more details see also: Henyey (1964) ApJ 139, 309
c-----------------------------------------------------------------------
c          This code is a slight adaption of Tjeerd's 1980 version:
c             - no plotting (subroutine PLORB + various komplot calls)
c             - some conversions to ANSI Fortran-77 + some '90' extensions:
c               (generic names, such as SQRT, LOG etc.)
c               (IF/THEN/ELSE/ENDIF, DOWHILE/ENDDO, DO/ENDDO)
c             - common blocks in include files (henyey.h, maxorb.h, model.h)
c               a parameter MAXORB is used for array sizes
c             - disk parameters are also in the /MODEL/ common
c             - each routine has some kind of IMPLICIT REAL***
c             - ``variable'' declaration of "real*8" which can al be
c               easily turned into "real*4": (note use of lower case here!)
c                   sed 's/REAL\*8/REAL\*4/gp' $file
c               or:
c                   sed 's/REAL\*4/REAL\*8/gp' $file
c		(note use of upper case here!)
c               The VIPDA routine usesa special declaration
c               DOUBLE PRECISION, which should remain that way!
c		There are however also locations in the code where
c		manual changes are needed ("metot.LE.1e-" and
c		"DATA eps/1e-")
c		==> Those locations are marked with the comment "FORCON"
c             - interface to NEMO's GETPARAM user-interface (ftoc - see
c               local Makefile how to do this)
c
c-----------------------------------------------------------------------
c Calling Structure of the program:
c
c MAIN--userinp
c       |
c       start---diffun--bar1
c       |       |       |
c       |       |       disk2
c       |       |
c       |       difs----diffun--bar1
c       |                       |
c       |                       disk2
c       |
c       orbit---henyey--cstruc--bar1
c                       |       |
c                       |       disk2
c                       |
c                       glesom--decom---vipda
c                               |
c                               fbsumb--vipda
c-----------------------------------------------------------------------
c  History:
c     17-sep-91   created - and sort of works         Peter Teuben (PJT)
c     18-sep-91   NEMO user interface - some workings unclear       PJT
c     19-sep-91   fixed up a few more old style fortrans            PJT
c     24-sep-91   toyed with the newly checkly C version of BARX    PJT
c                 dump (I,t/T,X,Y) of orbits on LUN=20 - for mongo
c     23-dec-23   trying to cleanup                                 PJT
c     
c***********************************************************************
c Delete comments 'cm' and fiddle with USERINP to make it standalone
cm      PROGRAM main
cm      CALL nemomain
cm      END
c***********************************************************************
      SUBROUTINE nemomain
      IMPLICIT REAL*8(a-z)
      include 'henyey.h'
      INTEGER n,n1,n2,l
      LOGICAL out
c                   (LUN1 is for input, LUN2 for output) 
      lun1=5
      lun2=6
c                   (Get stuff from the commandline via NEMO)
      CALL userinp(1)
c                   (some constants, kept in the /CONSTS/ common in model.h)
      pi=3.14159265358979
      gravc=0.4298
      gravc=1
      ome2=ome*ome
c                   (MASS, A, B are those of the bar)
      mass = 2.29
      mass=mass*1.33333333/2.0
      a=0.5
      b=0.01
      
c                   MDi mass of a disk, ADi length scale of a disk (i=1,2)
c                   Note MDi is actually M*GRAVC
      md1=9.9
      ad1=1.414
      md2=1.23
      ad2=0.106
c--BEGIN TEST - reset some variables for known analytical answers
      IF(.TRUE.) THEN
        mass = 1.0
        a = 1.0
        b = 0.2
        md1=0.0
        ad1=2.0
        md2=0.0
        ad2=1.0
      ENDIF
c--END TEST
c
      n1=MAXORB-1
      n2=n1+1
c
      DO n=1,n2
         tau(n)=FLOAT(n-1)/FLOAT(n1)
      ENDDO
c
c        Estimate values for parameters of periodic orbit
c        with a trial integration
c
      CALL start(n1)
c
      out=.TRUE.            
      uit=.FALSE.
      CALL orbit(n1,out)
c                           Find more orbits ...(these are for type 1)
      DO l=1,norbits
         IF(type.EQ.1) THEN
            x0=x0+step
            y0=0.0
            u0=0.0
            v0=v(1)
         ELSE IF (type.EQ.2.OR.type.EQ.3) THEN
            x0=0.0
            y0=y0+step
            u0=u(1)
            v0=0.0
         ELSE
            WRITE(lun2,'(''Illegal type'')')            
            STOP
         ENDIF
         CALL orbit(n1,out)
      ENDDO
c
      END
c***********************************************************************      
      SUBROUTINE orbit(n1,out)
      IMPLICIT REAL*8(a-z)
      INCLUDE 'henyey.h'
      INTEGER n,n1,n2,tel,i
      LOGICAL out,dum
      REAL*8 xi(MAXORB),yi(4,MAXORB),ea(4,MAXORB)
      REAL*8 zk(2),ap(50),me(4),syi(4)
c
      WRITE(lun2,2001)
      WRITE(lun2,2005) x0,y0,u0,v0,type,ome
      WRITE(lun2,2001) 
      
      dum=out
      tel=0
      n2=n1+1
   10 CALL henyey(4,1,n1,xi,yi,zk,ap,ea,1)
      tel=tel+1
      IF(tel.EQ.60) GOTO 200
      DO i=1,4
         me(i)=0.0
         syi(i)=0.0
         DO n=1,n2
            me(i) = me(i) + ABS(ea(i,n))
            syi(i)= syi(i)+ ABS(yi(i,n))
         ENDDO
         me(i)=me(i)/syi(i)
      ENDDO
      metot=me(1)+me(2)+me(3)+me(4)
      IF(metot-0.002) 21,21,20
   20 ucy=0.5
      IF(metot.GT.0.4) ucy=0.1
      GOTO 30
   21 ucy=1.0
   30 CONTINUE
      per=per+ucy*(zk(1)-per)
      WRITE(lun2,2010) me(1),me(2),me(3),me(4),ucy,per,u(1),v(1)
      DO n=1,n2
         x(n)=x(n)+ucy*ea(1,n)
         y(n)=y(n)+ucy*ea(2,n)
         u(n)=u(n)+ucy*ea(3,n)
         v(n)=v(n)+ucy*ea(4,n)
      ENDDO
c     FORCON Note: for double precision you may need smaller test value
cREAL4:
c      IF(metot.LT.1e-6) GOTO 100
cREAL8:
      IF(metot.LT.1d-9) GOTO 100
      GOTO 10
  100 IF(.NOT.out) GOTO 200
      uit=out
      out=.FALSE.
      GOTO 10
  200 out=dum
      uit=.FALSE.
      WRITE(lun2,2000)
c
c 2000:  newpage
c 2001:  newline

c	2000: newpage
c	2001: newline
 2000 FORMAT(" ")
 2001 FORMAT(1x)
 2005 FORMAT("x0",f12.8,8x,"y0",f12.8,8x,"u0",f12.8,8x,"v0",f12.8,
     *       8x,"type",i8,4x,"ome",f6.2)
 2010 FORMAT("dev",2x,4e10.3,4x,"ucy",f4.1,4x,"per",f16.10,4x,
     *       "u(1)",f16.10,4x,"v(1)",f16.10)
      
      END
c***********************************************************************      
      SUBROUTINE henyey(i1,k1,n1,xi,yi,zk,ap,ea,proc)
      IMPLICIT REAL*8(a-h,o-z)
      INTEGER proc
      include 'maxorb.h'
      INTEGER ipr(10)
      REAL*8 xi(MAXORB),yi(4,MAXORB),ea(4,MAXORB)
      REAL*8 fn(4,MAXORB),fd(4,7,MAXORB)
      REAL*8 zk(2),ap(50),ag(6,7),bg(6,7),dg(6),aa(7),bb(6,7),rm(6,7)
c
      n2=n1+1
      i2=i1+k1
      m1=i2+1
c
      DO j=1,m1
         DO i=1,i2
            bg(i,j)=0.0
            ag(i,j)=0.0
         ENDDO
         DO i=1,i1
            DO n=1,n2
               fd(i,j,n)=0.0
            ENDDO
         ENDDO
      ENDDO
c
      IF(proc.EQ.1) CALL cstruc(n1,xi,yi,zk,ap,fn,fd,ag,bg,dg)
c
      DO n=1,n1
         b=-0.5*(xi(n+1)-xi(n))
         DO i=1,i1
            DO j=1,i1
               bb(i,j)=-b*fd(i,j,n)
            ENDDO
            bb(i,i)=bb(i,i)+1.0
            DO j=1,i2
               rm(i,j)=b*(fd(i,j,n+1)+fd(i,j,n))
            ENDDO
            rm(i,m1)=yi(i,n+1)-yi(i,n)+b*(fn(i,n+1)+fn(i,n))
         ENDDO
         
         CALL glesom(bb,6,i1,ipr,rm,m1,err)
         
         DO i=1,i1
            DO j=1,m1
               fd(i,j,n)=rm(i,j)
            ENDDO
         ENDDO
      ENDDO
c
c  Recurrence relations
c
      DO i=1,i2
         ag(i,m1)=dg(i)
         DO 41 m=1,i1
            IF(ag(i,m))42,41,42
   41    CONTINUE
         GOTO 45
   42    DO n=1,n1
            DO l=1,i1
               aa(l)=ag(i,l)
            ENDDO
            DO j=1,m1
               DO l=1,i1
                  ag(i,j)=ag(i,j)+aa(l)*fd(l,j,n)  
               ENDDO
            ENDDO
         ENDDO
   45    CONTINUE   
         DO j=1,i1
               ag(i,j)=ag(i,j)+bg(i,j)
         ENDDO
         rm(i,1)=-ag(i,m1)
      ENDDO
c
      CALL glesom(ag,6,i2,ipr,rm,1,err)
c
c Improved solution
c
      aa(m1)=1.0
      m=n2
      DOWHILE(m.GT.0)
         DO i=1,i1
            ea(i,m)=rm(i,1)
         ENDDO
         IF(m.GT.1) THEN
            DO i=1,i2
               aa(i)=rm(i,1)
            ENDDO
            DO i=1,i1
               DO j=1,m1
                  rm(i,1)=rm(i,1)+fd(i,j,m-1)*aa(j)
               ENDDO
            ENDDO
         ENDIF
         m=m-1
      ENDDO

      IF (k1.GT.0) THEN
         DO k=1,k1
            l=i1+k
            zk(k)=zk(k)+rm(l,1)
         ENDDO
      ENDIF
      END
c***********************************************************************      
      SUBROUTINE glesom(a,nr,n,ipr,b,m,d1)
      IMPLICIT REAL*8(a-h,o-z)
c
c  To solve Ax=b having M right-hand sides by LU decomposition (DECOM)
c  and forward and backward substitution (FBSUBM)
c
      REAL*8 a(nr,1),b(nr,1),work(10)
      INTEGER ipr(1)
      EXTERNAL vipda
c      
c                    The decomposition A=LU
c===likely BUG: 4th argument of DECOM is wrong, must be a real 1D array? - PJT
      CALL decom(a,nr,n,ipr,ipr,d1,vipda)     
c      CALL decom(a,nr,n,work,ipr,d1,vipda)     
c
c                    Tests if A appears singular, if not
c                    solve AX=B by forward and backward substitutions
c
      IF(d1.NE.0.0) THEN
         CALL fbsubm(a,nr,n,ipr,b,m,vipda)
      ENDIF
      END
c***********************************************************************      
      SUBROUTINE decom(a,nr,n,v,ipr,d1,vipda)
      IMPLICIT REAL*8(a-h,o-z)
c to decompose A square matrix into lower traingular and upper
c triangular
      REAL*8 A(nr,1),v(1)
      INTEGER ipr(1)
      EXTERNAL vipda
c     FORCON Note: for double precision you may want smaller eps
cREAL4:
c      DATA eps/1.0e-6/
cREAL8:
      DATA eps/1.0d-12/
      e8=8.0*eps
c                       calculate the Euclidean norm
      DO i=1,n
         y=0
         CALL vipda(a(i,1),nr,a(i,1),nr,n,y)
         IF(y.LE.0) THEN
            d1=0.0
            RETURN
         ENDIF
         v(i)=1.0/SQRT(y)
      ENDDO
      d1=1.0
      DO k=1,n
         l=k
         x=0.0
         k1=k-1
         DO i=k,n
            y=-a(i,k)
            IF(k1.GT.0) CALL vipda(a(i,1),nr,a(1,k),1,k1,y)
            a(i,k)=-y
            y=ABS(y*v(i))
            IF(x.LT.y) THEN            
               x=y
               l=i
            ENDIF
         ENDDO
         IF(l.NE.k) THEN
            d1=-d1
c                                Interchange the rows k and l
            DO j=1,n
               y=a(k,j)
               a(k,j)=a(l,j)
               a(l,j)=y
            ENDDO
            v(l)=v(k)
         ENDIF
         ipr(k)=l
c                                Test if appears singular
         IF(x.LT.e8) THEN
            d1=0.0
            RETURN
         ENDIF
c                                Compute elements of upper triangular 
c                                matrix with diagonal elements excluded
         x=-1.0/a(k,k)
         j=k+1
         DOWHILE(j.LE.n)
            y=-a(k,j)
            IF(k1.GT.0) CALL vipda(a(k,1),nr,a(1,j),1,k1,y)
            a(k,j)=x*y
            j=j+1
         ENDDO
      ENDDO
      END
c***********************************************************************            
      SUBROUTINE fbsubm(a,nr,n,ipr,b,m,vipda)
      IMPLICIT REAL*8(a-h,o-z)
c
c  To solve Ax=b having M right-hand sides using the results from
c  DECOM by forward and backward substitutions
      REAL*8 a(nr,1),b(nr,1)
      INTEGER ipr(1)
      EXTERNAL vipda

c                 Interchange elements of B
      DO i=1,n
         IF(ipr(i).NE.i) THEN
            DO k=1,m
               x=b(i,k)
               j=ipr(i)
               b(i,k)=b(j,k)
               b(j,k)=x
            ENDDO
         ENDIF
      ENDDO
c                 Solves LY=B by forward substitution
      DO k=1,m
         b(1,k)=-b(1,k)/a(1,1)
         IF(n.NE.1) THEN
            DO i=2,n
               x=b(i,k)
               CALL vipda(a(i,1),nr,b(1,k),1,i-1,x)
               b(i,k)=-x/a(i,i)
            ENDDO
         ENDIF
c                 Solves UX=Y by backward substitution
         i=n
         DOWHILE(i.GE.1)
            x=b(i,k)
            ni=n-1
            IF(ni.GT.0) THEN
               i1=i+1
               CALL vipda(a(i,i1),nr,b(i1,k),1,ni,x)
            ENDIF
            b(i,k)=-x
            i=i-1
         ENDDO
      ENDDO

      END
c***********************************************************************      
      SUBROUTINE vipda(a,na,b,nb,n,c)
      IMPLICIT REAL*8(a-h,o-z)
c
c  Inner product subroutine:
c    Array A is taken with a stride of NA
c    Array B is taken with a stride of NB
c  The product is done with length N.
c  Incoming value of C is added to the inner product and result
c  is stored in C also.
c
      REAL*8 a(1),b(1)
      DOUBLE PRECISION sum, da, db

      sum=c
      nn=1-nb
      jl=na*n
      DO jj=1,jl,na
         nn=nn+nb
         da=a(jj)
         db=b(nn)
         sum=sum+da*db
      ENDDO      
      c=sum
      END
c***********************************************************************      
      SUBROUTINE cstruc(n1,xi,yi,zk,ap,fn,fd,ag,bg,dg)
      IMPLICIT REAL*8(a-z)
      INCLUDE 'henyey.h'
      INTEGER n,n1,n2,i,j
      REAL*8 zk(2),ap(50),ag(6,7),bg(6,7),dg(6)
      REAL*8 xi(MAXORB),yi(4,MAXORB),fn(4,MAXORB),fd(4,7,MAXORB)
c
      n2=n1+1
      zk(1)=per
c                    Coefficients in linearized constraints
      DO i=1,6
         DO j=1,7
            ag(i,j)=0.0
            bg(i,j)=0.0
         ENDDO
      ENDDO
      ag(1,1)=1.0
      ag(2,2)=1.0
      dg(6)=0.0
c                    Handle different types of orbits      

c                        Type 1 orbit  (start on X axis, end on Y axis)
c                        Handles X_1, X_2, X_3 and X_4 orbits
      IF(type.EQ.1) THEN    
         ag(3,3)=1.0
         bg(4,1)=1.0
         bg(5,4)=1.0
         dg(1)=x(1)-x0
         dg(2)=y(1)
         dg(3)=u(1)
         dg(4)=x(n2)
         dg(5)=v(n2)
      ELSE IF (type.EQ.2) THEN
c                         Type 2 orbit (start on Y axis, end on X axis)
c                         Handles X_1, X_2, X_3 and X_4 orbits
         ag(3,4)=1.0
         bg(4,2)=1.0
         bg(5,3)=1.0
         dg(1)=x(1)
         dg(2)=y(1)-y0
         dg(3)=v(1)
         dg(4)=y(n2)
         dg(5)=u(n2)
      ELSE IF (type.EQ.3) THEN
c                         Type 3 orbit (start on Y axis, end on Y axis)
c                         Handles X_1, X_2, X_3 and X_4 orbits, but also
c                         the SPO and LPO around stable Lagrangian points
         ag(3,4)=1.0
         bg(4,1)=1.0
         bg(5,4)=1.0
         dg(1)=x(1)
         dg(2)=y(1)-y0
         dg(3)=v(1)
         dg(4)=x(n2)
         dg(5)=v(n2)
      ELSE
         WRITE(lun2,'(''Illegal type '')')
         STOP
      ENDIF
      CALL bar1(x0,y0,pgb,pxb,pyb,pxxb,pxyb,pyyb)
      CALL disk2(md1,ad1,x0,y0,pgd1,pxd1,pyd1,pxxd1,pxyd1,pyyd1)
      CALL disk2(md2,ad2,x0,y0,pgd2,pxd2,pyd2,pxxd2,pxyd2,pyyd2)
      pgd=pgd1+pgd2
      phi=pgb+pgd+0.5*ome2*(x0*x0+y0*y0)
      jaci1=0.5*(u(1)*u(1)+v(1)*v(1))-phi
c
      IF(uit) THEN
         WRITE(lun2,2001)
         WRITE(lun2,2005)
         WRITE(lun2,2001)
      ENDIF
      DO n=1,n2
         xx=x(n)
         yy=y(n)
         uu=u(n)
         vv=v(n)
         xi(n)=tau(n)
         yi(1,n)=xx
         yi(2,n)=yy
         yi(3,n)=uu
         yi(4,n)=vv
         CALL bar1(xx,yy,pgb,pxb,pyb,pxxb,pxyb,pyyb)
         CALL disk2(md1,ad1,xx,yy,pgd1,pxd1,pyd1,pxxd1,pxyd1,pyyd1)
         CALL disk2(md2,ad2,xx,yy,pgd2,pxd2,pyd2,pxxd2,pxyd2,pyyd2)
         pgd=pgd1+pgd2
         pxd=pxd1+pxd2         
         pyd=pyd1+pyd2
         pxxd=pxxd1+pxxd2
         pxyd=pxyd1+pxyd2
         pyyd=pyyd1+pyyd2
         phi=pgb+pgd+0.5*ome2*(xx*xx+yy*yy)
         phix=pxb+pxd+ome2*xx
         phiy=pyb+pyd+ome2*yy
         phixx=pxxb+pxxd+ome2
         phixy=pxyb+pxyd
         phiyy=pyyb+pyyd+ome2
         fn(1,n)=per*uu
         fn(2,n)=per*vv
         fn(3,n)=per*( 2.0*ome*vv+phix)
         fn(4,n)=per*(-2.0*ome*uu+phiy)
         IF(uit .AND. MOD(n,10).EQ.1) THEN
            rr=SQRT(xx*xx+yy*yy)
            cph=xx/rr
            sph=yy/rr
            vr=vv*sph+uu*cph
            vt=vv*cph-uu*sph
            uu0=uu-ome*yy
            vv0=vv+ome*xx
            vr0=vv0*sph+uu0*cph
            vt0=vv0*cph-uu0*sph
            vel2=uu*uu+vv*vv
            stab=(-phixx*vv*vv+phixy*2.0*uu*vv-phiyy*uu*uu)/vel2
     *           +3.0*(-2.0*ome+(uu*phiy-vv*phix)/vel2)**2
            ekin=0.5*vel2
            jaci=ekin-phi
            ree=(jaci-jaci1)/jaci
            WRITE(lun2,2010) n,tau(n),xx,yy,uu,vv,vr,vt,uu0,vv0,vr0,vt0,
     *                       phi,jaci,ree,stab
         ENDIF
         IF (uit) THEN
            WRITE(20,*) n,tau(n),xx,yy
         ENDIF
         DO i=1,4
            fd(i,6,n)=0.0
         ENDDO
         fd(1,5,n)=fn(1,n)/per
         fd(2,5,n)=fn(2,n)/per
         fd(3,5,n)=fn(3,n)/per
         fd(4,5,n)=fn(4,n)/per
         fd(1,1,n)=0.0
         fd(2,1,n)=0.0
         fd(3,1,n)=per*phixx
         fd(4,1,n)=per*phixy
         fd(1,2,n)=0.0    
         fd(2,2,n)=0.0             
         fd(3,2,n)=per*phixy
         fd(4,2,n)=per*phiyy
         fd(1,3,n)=per
         fd(2,3,n)=0.0
         fd(3,3,n)=0.0
         fd(4,3,n)=-2.0*ome*per
         fd(1,4,n)=0.0
         fd(2,4,n)=per
         fd(3,4,n)=2.0*ome*per
         fd(4,4,n)=0.0
      ENDDO
      DO i=1,50
         ap(i)=0.0
      ENDDO
      IF(uit) WRITE(lun2,2001)

c 2000: newpage
c 2001: newline
 2000 FORMAT(" ")
 2001 FORMAT(1x)
 2005 FORMAT(" n",3x,"tau",6x,"x",7x,"y",7x,"u",7x,"v",7x,"vr",6x,
     *       "vt",6x,"u0",6x,"v0",5x,"vr0",5x,"vt0",5x,
     *       "phi",5x,"jacobi",5x,"rel err",6x,"stab")
 2010 FORMAT(i3,f6.2,11f8.4,f12.7,2e10.2)
      END
c***********************************************************************            
      SUBROUTINE start(n1)
      IMPLICIT REAL*8(a-z)
      include 'henyey.h'
      INTEGER n,n1,n2,l
      REAL*8 yy(4),dy(4),yypr(4)
      REAL*8 t(MAXORB)
c                             Input of the first orbit: see now USERINP
      n2=n1+1
c                                
      WRITE(lun2,2005) x0,y0,u0,v0,per,type
 2005 FORMAT("x0",f12.5,8x,"y0",f12.8,8x,"u0",f12.8,8x,"v0",f12.8,
     *       8x,"per",f12.8,8x,"type",i8)
c  
      DO n=1,n2
         t(n)=FLOAT(n-1)*ABS(per)/FLOAT(n1)
      ENDDO     
c
      x(1)=x0
      y(1)=y0
      u(1)=u0
      v(1)=v0
      yy(1)=x(1)
      yy(2)=y(1)
      yy(3)=u(1)
      yy(4)=v(1)
      h=ABS(per)*0.0005*0.9999
c
      tt=0.0
      CALL diffun(yy,dy)

      n=2
      DOWHILE(n.LE.n2)
         DOWHILE(tt.LT.t(n))
            tprev=tt
            DO l=1,4
               yypr(l)=yy(l)
            ENDDO
            CALL difs(tt,yy,dy,h)
         ENDDO
         fr=(t(n)-tprev)/(tt-tprev)
         x(n)=yypr(1)+fr*(yy(1)-yypr(1))
         y(n)=yypr(2)+fr*(yy(2)-yypr(2))
         u(n)=yypr(3)+fr*(yy(3)-yypr(3))
         v(n)=yypr(4)+fr*(yy(4)-yypr(4))
         n=n+1
      ENDDO

      END
c***********************************************************************            
      SUBROUTINE difs(t,y,dy,h)
      IMPLICIT REAL*8(a-h,o-z)
      REAL*8 y(4),dy(4)
      
      CALL diffun(y,dy)
      DO i=1,4
         y(i)=y(i)+h*dy(i)
      ENDDO
      t=t+h
      
      END
c***********************************************************************            
      SUBROUTINE diffun(yy,dyn)
      IMPLICIT REAL*8(a-z)
      REAL*8 yy(4),dyn(4)
      include 'model.h'

      CALL bar1(         yy(1),yy(2),pgb, pxb, pyb, pxxb, pxyb, pyyb)
      CALL disk2(md1,ad1,yy(1),yy(2),pgd1,pxd1,pyd1,pxxd1,pxyd1,pyyd1)
      CALL disk2(md2,ad2,yy(1),yy(2),pgd2,pxd2,pyd2,pxxd2,pxyd2,pyyd2)
      pxd=pxd1+pxd2
      pyd=pyd1+pyd2
      phix=pxb+pxd+ome2*yy(1)
      phiy=pyb+pyd+ome2*yy(2)
      dyn(1)=yy(3)
      dyn(2)=yy(4)
      dyn(3)= 2.0*ome*yy(4)+phix
      dyn(4)=-2.0*ome*yy(3)+phiy
      END
c***********************************************************************       
      SUBROUTINE bar1(x,y,pot,fx,fy,fxx,fxy,fyy)
      IMPLICIT REAL*8(a-z)
      include 'model.h'
c
c       Homogeneous prolate spheroid with axes A,B and mass MASS
c       Potential is defined > 0
c       BUG: can't handle spherical (a.EQ.b) cases
c
c                           If zero mass, quick exit
      IF(mass.EQ.0.0) THEN
         pot=0.0
         fx=0.0
         fy=0.0
         fxx=0.0
         fxy=0.0
         fyy=0.0
         RETURN
      ENDIF
      gm=mass*gravc
      IF (.TRUE.) THEN
c				testing a C routine BARX (see any barX.c)
         CALL barx(gm,a,b,x,y,pot,fx,fy,fxx,fxy,fyy)
         RETURN
      ENDIF
      ecc=SQRT(a*a-b*b)/a
      ae=a*ecc
      ae3=ae*ae*ae
      ale=LOG((1.0-ecc)/(1.0+ecc))
      aint1=(-2.0*ecc-ale)/ae3
c      omcr=SQRT(1.5*gm*aint1)         ! redundant for now ...
      r1=SQRT((x+ae)*(x+ae)+y*y)      
      r2=SQRT((x-ae)*(x-ae)+y*y)      
      IF(r1+r2.LE.2.0*a) THEN
c                                         (x,y) interior point
         al=0.0
         ale=LOG((1.0-ecc)/(1.0+ecc))
         aint1=(-2.0*ecc-ale)/ae3
         aint2=(ecc/(1.0-ecc*ecc)+0.5*ale)/ae3
         aint3=-ale/ae
         aux1=0.0
         aux2=0.0
         aux3=0.0
c                                         (x,y) exterior point
      ELSE
         bb=a*a+b*b-x*x-y*y
         cc=a*a*b*b-b*b*x*x-a*a*y*y
         al=-0.5*bb+0.5*SQRT(bb*bb-4.0*cc)
         rtal=SQRT(a*a+al)
         ale=LOG((rtal-ae)/(rtal+ae))
         aint1=(-2.0*ae/rtal-ale)/ae3
         aint2=(ae*rtal/(b*b+al)+0.5*ale)/ae3
         aint3=-ale/ae
         daldx=2.0*x*(al+b*b)/(2.0*al+bb)
         daldy=2.0*y*(al+a*a)/(2.0*al+bb)
         aux1=1.5*gm*x*daldx/(rtal*(a*a+al)*(b*b+al))
         aux2=1.5*gm*y*daldy/(rtal*(b*b+al)*(b*b+al))
         aux3=1.5*gm*x*daldy/(rtal*(a*a+al)*(b*b+al))
      ENDIF
      pot=-x*x*aint1-y*y*aint2+aint3
      pot=0.75*gm*pot
      fx=-1.5*gm*x*aint1
      fy=-1.5*gm*y*aint2
      fxx=-1.5*gm*aint1+aux1
      fyy=-1.5*gm*aint2+aux2
      fxy=              aux3
      
      END
c***********************************************************************           
      SUBROUTINE disk2(gm,a,x,y,psi,psix,psiy,psixx,psixy,psiyy)
      IMPLICIT REAL*8(a-z)
c
c Toomre disk no. 2
c gm = gravc*mass     a = scale length
c potential psi > 0   fx = dspidx = psix
c   BUG: can't handle the origin
c
c                           If zero mass, quick exit
      IF(gm.EQ.0.0) THEN
         psi=0.0
         psix=0.0
         psiy=0.0
         psixx=0.0
         psixy=0.0
         psiyy=0.0
      ENDIF
      r2=x*x+y*y
      r=SQRT(r2)
      r3=r*r2
      a2=a*a
      a2r2=a2+r2
      wa2r2=SQRT(a2r2)
      psi=gm/wa2r2
      psir=-gm*r/(a2r2*wa2r2)
      psirr=-gm*(a2-2.0*r2)/(a2r2*a2r2*wa2r2)
      psix=psir*x/r
      psiy=psir*y/r
      psixx= psir*y*y/r3+psirr*x*x/r2
      psixy=-psir*x*y/r3+psirr*x*y/r2
      psiyy= psir*x*x/r3+psirr*y*y/r2
      END
c***********************************************************************            
      SUBROUTINE userinp(n)
      IMPLICIT REAL*8(a-h,o-z)
      INTEGER n
c
c  This routine obtains parameters from the command line and
c  places them in the common blocks for later access 
c
      include 'henyey.h'
      DOUBLE PRECISION getdparam
      INTEGER getiparam

      WRITE(lun2,'(''REAL*8 version'')')
      IF (n.EQ.1) THEN
c                           NEMO user interface (see also 'getparam.inc')
         x0 = getdparam('x0')
         y0 = getdparam('y0')
         u0 = getdparam('u0')
         v0 = getdparam('v0')
         per = getdparam('per')
         type = getiparam('type')
         ome = getdparam('ome')
         norbits = getiparam('norbit')
         step = getdparam('step')
      ELSE
c                           do your own READ and WRITE here...
         WRITE(lun2,'(''Enter: x0,y0,u0,v0,per,type,ome,norbit,step'')')
         READ(lun1,*) x0,y0,u0,v0,per,type,ome,norbits,step
      ENDIF

      END
