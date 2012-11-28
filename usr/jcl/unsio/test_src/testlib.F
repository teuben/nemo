      program testlib
      implicit none

      integer lnblnk, narg, uns_init, uns_load, valid, ident,
     $     nbody, status, uns_sim_type, uns_get_eps
      character arg1 * 180, arg2 * 180, arg3 * 180
      real *4 eps

      narg = command_argument_count()
      if (narg.ne.3) then
         write(0,*) "You must give 3 parameters:"
         write(0,*) "filename selected_component selected_time"
         stop
      endif
      ! format fortran string to c string
      call get_command_argument(1,arg1)
      ! fortran to C string
      arg1(lnblnk(arg1)+1:lnblnk(arg1)+1) = '\0'
      print *, "filename :",arg1

      call get_command_argument(2,arg2)
      ! fortran to C string
      arg2(lnblnk(arg2)+1:lnblnk(arg2)+1) = '\0'

      call get_command_argument(3,arg3)
      ! fortran to C string
      arg3(lnblnk(arg3)+1:lnblnk(arg3)+1) = '\0'

      ! initialyze UNS engine
      ident=uns_init(arg1, arg2, arg3) ! return identifier for the snaphot
                                       ! ident must be positive
      write(0,*) "Uns type = ",uns_sim_type(ident)
      write (0,*) "UNS_INIT return :",ident
      if (ident.gt.0) then
         valid = 1
         do while (valid .gt. 0) ! loop on all the time step
            valid = uns_load(ident) ! load data belonging to ident snapshot
            write(0,*) "uns_load return=>",valid
            if (valid .gt. 0) then  ! it goes fine 
               call uns_get_nbody(ident, nbody) ! get #bodies
               write (0,*) ">>>> nbody =", nbody
               call start(ident,nbody)
               status=uns_get_eps(ident,"gas\0",eps);
               write(0,*) "gas eps=",eps
               status=uns_get_eps(ident,"halo\0",eps);
               write(0,*) "halo eps=",eps

            endif
         enddo
      endif
      end
      

      subroutine start(ident,nbody)
      implicit none
! input parameters
      integer ident,nbody
! UNS variable
      integer status
      integer  uns_get_range, uns_get_age, uns_get_metal, uns_get_cod
      real * 4 uns_get_eps, eps, rs
      real *4 time, pos(3,nbody), vel(3,nbody), mass(nbody), age(nbody),
     $     metal(nbody), tcod(7)
! various      
      integer io_nemo_f
      character * 90 out
      integer n,first,last,i
      write(0,*) '-1 '
      call uns_get_time(ident,time      )
      write (0,*) '>>>>>>> time =',time
      call uns_get_redshift(ident,rs      )
      write (0,*) '>>>>>>> RedShift =',rs
      
      call uns_get_pos (ident,pos ,nbody)
      call uns_get_vel (ident,vel ,nbody)
      call uns_get_mass(ident,mass,nbody)
      i=uns_get_cod(ident,'gas\0',time,tcod)
      if (i.eq.1) then 
         write (0,*) '>>>>>> cod gas i=',i
         write (0,*) (tcod(i), i=1,7)
      endif
      i=uns_get_cod(ident,'all\0',time,tcod)
      if (i.eq.1) then 
         write (0,*) '>>>>>> cod all i=',i
         write (0,*) (tcod(i), i=1,7)
      endif
      write(0,*) '-2 '
      status = uns_get_range(ident,"disk\0",n,first,last);
      write(0,*) '-22 '
      status = uns_get_range(ident,"gas\0",n,first,last);
      write(0,*) '-23 '
      if (status.gt.0) then
         n = uns_get_metal(ident,metal,nbody)
         write(0,*) 'Got',n,' particles with metalicity'
      endif
      write(0,*) '-24 '
      status = uns_get_range(ident,"stars\0",n,first,last);
      if (status.gt.0) then
         n = uns_get_age(ident,age,nbody)
         write(0,*) 'Got',n,' particles with ages'
      endif
      write(0,*) '-3 '
      status = uns_get_range(ident,"halo\0",n,first,last);
      status = uns_get_range(ident,"bulge\0",n,first,last);
      write(0,*) '-4 '


      out = "nemo.out"
      write(0,*) "nbody=",nbody," time=",time
      status = io_nemo_f(out,80,nbody,"float,save,3n,n,t,x,v,m\0",
     $     nbody,time,pos,vel,mass)
      
      end








