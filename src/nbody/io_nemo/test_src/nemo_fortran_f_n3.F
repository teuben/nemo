      program nemo_fortran
      implicit none

c maximum bodies
/*#define MAXBODY 100000*/
      integer MAXBODY
      parameter (MAXBODY=100000)

      real*4
c particles masses.
     +     mass(1:MAXBODY),
c particles positions and velocities    
     +     pos(MAXBODY,3), vel(MAXBODY,3),
c particles potential and accelerations
     +     pot(MAXBODY), acc(MAXBODY,3),
c snapshot time
     +     ts

      integer
     $     keys(MAXBODY)

c files name
      character infile*80,outfile*80, st*80, sp*80

      integer i,j,k,nbody,io_nemo_f,close_io_nemo_f

c get input snapshot filename
      write(*,*)'Infile name : '
      read(*,'(a80)')infile

c get output snapshot filename
      write(*,*)'outfile name : '
      read(*,'(a80)')outfile

c get selected time
      write(*,*)'selected time (put a # at the end) : '
      read(*,'(a80)')st

c get selected particles
      write(*,*)'selected particles (put a # at the end) : '
      read(*,'(a80)')sp

      i = 1
      do while (i.ne.0) 

c read the snapshot until there is not time steps ( i <= 0 )
         i=io_nemo_f(infile,80,MAXBODY,"real4,n3,read,n,m,x,v,t,k,
     $        info,p,st,sp,acc\0",
     $        nbody,mass,pos,vel,ts,keys,pot,st,sp,acc)


c save the snapshot 
         if (i.gt.0) then
            j=io_nemo_f(outfile,80,MAXBODY,"real4,n3,save,n,m,x,v,t,k,
     $           info,p,acc\0",nbody,mass,pos,vel,ts,keys,pot,acc)


         endif

      end do

      k= close_io_nemo_f(outfile,80)


      end
