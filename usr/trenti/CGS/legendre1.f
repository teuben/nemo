c$$$c
c$$$c
c$$$c
c$$$c  test driver
c$$$c
c$$$      
c$$$      program test
c$$$
c$$$      implicit none
c$$$      include 'common.blk'        !contains lmax
c$$$      double precision plm(lmax,lmax) !matrix with p_l^m(cos(theta)) 
c$$$      double precision dplm(lmax,lmax) !matrix with d/dx(p_l^m(x))
c$$$      common /pol/ plm,dplm     !please init all to zero
c$$$                                !note that P00 is stored in plm(1,1),
c$$$                                !P10 in plm(2,1)... and so on
c$$$
c$$$      double precision x,aux, plgndr
c$$$      external  plgndr
c$$$      integer l,i,j
c$$$      
c$$$      l=6
c$$$      x=.999d0
c$$$
c$$$      do i=1,lmax
c$$$         do j=1,lmax
c$$$           plm(j,i)=0.
c$$$           dplm(j,i)=0.
c$$$         enddo
c$$$      enddo
c$$$
c$$$            
c$$$      call  plgndrMT(l,x)
c$$$
c$$$c writes output
c$$$       do i=1,lmax
c$$$          write(*,*) 'l = ', i-1
c$$$         do j=1,lmax
c$$$            if(j.le.i) then
c$$$               aux= plgndr(i-1,j-1,x)
c$$$            else
c$$$               aux=0.
c$$$            endif
c$$$            
c$$$            write(*,*) j-1,plm(i,j)-aux,dplm(i,j),aux
c$$$         enddo
c$$$         write(*,*) '*****************'
c$$$      enddo
c$$$      
c$$$      write(*,*) plm(1,1)
c$$$
c$$$
c$$$      end
c$$$










c
c
c  prototipo di funzione che calcola i polinomi di legendre 
c  e le loro derivate fino ad ordine l dato un x
c
c  sfrutta la relazione di ricorrenza dal paragrafo 6.8 di Numerical recipes
c
c


      subroutine plgndrMT(l,x)
      
      implicit none
      
      include 'common.blk'        !contains lmax

      integer l                 !spherical armonics truncation
      double precision x        !argument (x=cos(theta))
      
      double precision plm(lmax,lmax) !matrix with p_l^m(cos(theta)) 
      double precision dplm(lmax,lmax) !matrix with d/dx(p_l^m(x))
      common /pol/ plm,dplm     !please init all to zero
                                !note that P00 is stored in plm(1,1),
                                !P10 in plm(2,1)... and so on

c internal variables
      integer m,i,ll            !loop indexes
      double precision pmm,pmm1,pll,pmmp1 !aux leg polinomials
      double precision fact     !for odd factorial computation 
      double precision dpmm,dpmm1,dpll,dpmmp1 !aux derivatives leg polinomials
      double precision somx2,somx22 !arguments for starting values of Plm(x)
      

cccccccccccccccccccccccc
c legendre polinomials c
cccccccccccccccccccccccc


      
c     argument assignment
      somx22=(1.-x)*(1.+x)
      somx2=dsqrt(somx22)
      

c loop over m
      do m=0,l

         pmm=1.d0                 !compute P^M_M
         fact=1.
         do i=1,m
            pmm=-pmm*fact*somx2
            fact=fact+2.d0
         enddo
         plm(m+1,m+1)=pmm             !matrix assignment

         pmmp1=x*(2*m+1)*pmm      !compute P^M_(M+1)
         plm(m+2,m+1)=pmmp1           !matrix assignment
         
         do ll=m+2,l              !compute P^M_(L)
            pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
            pmm=pmmp1
            pmmp1=pll
            plm(ll+1,m+1)=pll         !matrix assignment
         enddo

      enddo

ccccccccccccccccccccccccccccccccccccccc
c derivatives of legendre polinomials c
ccccccccccccccccccccccccccccccccccccccc
      
      dplm(1,1)=0.
      do ll=1,l                 !analytical relation
         do m=0,ll
            dplm(ll+1,m+1)=-(-(m+ll)*plm(ll,m+1)+ll*x*plm(ll+1,m+1))
     &           /somx22  
         enddo
      enddo

      return
      end



C
C
C same as  plgndrMT(l,x) but without derivatives calculation
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine plgndrMT1(l,x)
      
      implicit none
      
      include 'common.blk'        !contains lmax

      integer l                 !spherical armonics truncation
      double precision x        !argument (x=cos(theta))
      
      double precision plm(lmax,lmax) !matrix with p_l^m(cos(theta)) 
      double precision dplm(lmax,lmax) !matrix with d/dx(p_l^m(x))
      common /pol/ plm,dplm     !please init all to zero
                                !note that P00 is stored in plm(1,1),
                                !P10 in plm(2,1)... and so on

c internal variables
      integer m,i,ll            !loop indexes
      double precision pmm,pmm1,pll,pmmp1 !aux leg polinomials
      double precision fact     !for odd factorial computation 
      double precision dpmm,dpmm1,dpll,dpmmp1 !aux derivatives leg polinomials
      double precision somx2,somx22 !arguments for starting values of Plm(x)
      

cccccccccccccccccccccccc
c legendre polinomials c
cccccccccccccccccccccccc

      
c     argument assignment
      somx22=(1.-x)*(1.+x)
      somx2=dsqrt(somx22)
      

c loop over m
      do m=0,l

         pmm=1.d0                 !compute P^M_M
         fact=1.
         do i=1,m
            pmm=-pmm*fact*somx2
            fact=fact+2.d0
         enddo
         plm(m+1,m+1)=pmm             !matrix assignment

         pmmp1=x*(2*m+1)*pmm      !compute P^M_(M+1)
         plm(m+2,m+1)=pmmp1           !matrix assignment
         
         do ll=m+2,l              !compute P^M_(L)
            pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
            pmm=pmmp1
            pmmp1=pll
            plm(ll+1,m+1)=pll         !matrix assignment
         enddo

      enddo

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      FUNCTION plgndr(l,m,x)
      INTEGER l,m
      DOUBLE PRECISION plgndr,x
      INTEGER i,ll
      DOUBLE PRECISION fact,pll,pmm,pmmp1,somx2
      if(m.lt.0.or.m.gt.l.or.abs(x).gt.1.)pause
     *'bad arguments in plgndr'
      pmm=1.
      if(m.gt.0) then
        somx2=sqrt((1.-x)*(1.+x))
        fact=1.
        do 11 i=1,m
          pmm=-pmm*fact*somx2
          fact=fact+2.
11      continue
      endif
      if(l.eq.m) then
        plgndr=pmm
      else
        pmmp1=x*(2*m+1)*pmm
        if(l.eq.m+1) then
          plgndr=pmmp1
        else
          do 12 ll=m+2,l
            pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
            pmm=pmmp1
            pmmp1=pll
12        continue
          plgndr=pll
        endif
      endif
      return
      END
