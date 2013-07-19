      program unsio_fortran
      implicit none

      integer lnblnk, narg, uns_init, uns_load, valid, ident,
     $     nbody,  uns_sim_type, cpt, uns_get_value_i,
     $     status
      !external iargc,getarg
      character ain * 80, aout * 80, atyp * 80, acomp * 80, atim * 80
      character out * 80, ii * 6, simdir * 100

      !narg = iargc()
      narg = command_argument_count()
      write(0,*) 'narg=',narg
      if (narg.ne.5) then
          call getarg(0,ain)
         write(0,*) "You must give 5 parameters:"
         write(0,*) ain(1:lnblnk(ain)),
     $        " inputfile outputfile selected_component",
     $        "  outttype selected_time"
          write(0,*) ""
          write(0,*) "Example : "

          write(0,*) ain(1:lnblnk(ain)),
     $         " sgs029 sgs029.gas gas gadget2 all"
         stop
      endif

      ! format fortran string to c string
      call get_command_argument(1,ain)

      call get_command_argument(2,aout)

      call get_command_argument(3,acomp)

      call get_command_argument(4,atyp)

      call get_command_argument(5,atim)

      ! initialyze UNS engine
      ident=uns_init(ain, acomp, atim) ! return identifier for the snaphot
                                       ! ident must be positive
      write(0,*) "Uns type = ",uns_sim_type(ident)
      write (0,*) "UNS_INIT return :",ident
      if (ident.gt.0) then
         valid = 1
         cpt = 0
         call uns_sim_dir(ident,simdir)
         write (0,*) "Simulation DIR = [",simdir(1:lnblnk(simdir)),"]"
         do while (valid .gt. 0) ! loop on all the time step
            valid = uns_load(ident) ! load data belonging to ident snapshot
            write(0,*) "uns_load return=>",valid
            if (valid .gt. 0) then  ! it goes fine 

               !! build output file name
               write (ii,'(I6)') cpt
               out = trim(aout)//"."//trim(adjustl(ii))
               write (0,*) out
               !! get nbody
               status=uns_get_value_i(ident,"nsel", nbody) ! get #bodies
               write (0,*) "nbody =", nbody
               !! start main subroutine
               call start(ident,nbody,out,atyp)
               cpt = cpt+1
            endif
         enddo
      endif
      end
      

      subroutine start(ident,nbody,out,atyp)
      implicit none
! input parameters
      integer ident,nbody
! UNS variable
      integer status
      integer  uns_get_range, uns_save_init
      integer  uns_get_value_f, uns_get_array_f, nsel
      integer  uns_set_array_f, uns_set_value_f, ok
      integer  uns_get_array_i
      real *4 time, pos(3,nbody), vel(3,nbody), mass(nbody), rho(nbody),
     $     hsml(nbody), u(nbody)
      integer id(nbody)
      character out * 80, atyp * 80
! various      
      integer n,first,last, idento

      ! initialyze uns output
      idento=uns_save_init(out,atyp)

      status=uns_get_value_f(ident,"time",time      )
      write (0,*) 'time =',time

      
      ! prepare time for saving
      ok=uns_set_value_f(idento,"time", time)

      ! load pos vel mass
      nsel=uns_get_array_f(ident,"all","pos" ,pos ,nbody)
      nsel=uns_get_array_f(ident,"all","vel" ,vel ,nbody)
      nsel=uns_get_array_f(ident,"all","mass",mass,nbody)
      write (0,*) 'nsel =',nsel,' nbody=',nbody
      nsel=uns_get_array_i(ident,"all","id",id,nbody)
      write (0,*) 'ID nsel =',nsel,' id(1)',id(1)
      ! check if gas component exist
      status = uns_get_range(ident,"bndry",n,first,last);

      status = uns_get_range(ident,"gas",n,first,last);
      if (status.gt.0) then ! prepare gas data for saving
         write (0,*) 'Gas ok, n=',n,' first=',first, ' last=',last
         ok=uns_set_array_f(idento,"gas","mass",mass(first)  ,n);
         ok=uns_set_array_f(idento,"gas","pos" ,pos (1,first),n);
         ok=uns_set_array_f(idento,"gas","vel" ,vel (1,first),n);

         ! gas density
         n=uns_get_array_f(ident,"gas","rho",rho,nbody);
         if (n.gt.0) then
            ok=uns_set_array_f(idento,"gas","rho" ,rho,n);
         endif

         ! gas hsml
         n=uns_get_array_f(ident,"gas","hsml",hsml,nbody);
         if (n.gt.0) then
            ok=uns_set_array_f(idento,"gas","hsml" ,hsml,n);
         endif                  
         
         ! gas internal energy
         n=uns_get_array_f(ident,"gas","u",u,nbody);
         if (n.gt.0) then
            ok=uns_set_array_f(idento,"gas","u" ,u,n);
         endif                  

      endif
      ! check if disk component exist
      status = uns_get_range(ident,"disk",n,first,last);
      if (status.gt.0) then ! prepare disk data for saving
         ok=uns_set_array_f(idento,"disk","mass",mass(first)  ,n);
         ok=uns_set_array_f(idento,"disk","pos" ,pos (1,first),n);
         ok=uns_set_array_f(idento,"disk","vel" ,vel (1,first),n);
      endif
      ! check if stars component exist
      status = uns_get_range(ident,"stars",n,first,last);
      if (status.gt.0) then ! prepare stars data for saving
         ok=uns_set_array_f(idento,"stars","mass",mass(first)  ,n);
         ok=uns_set_array_f(idento,"stars","pos" ,pos (1,first),n);
         ok=uns_set_array_f(idento,"stars","vel" ,vel (1,first),n);
      endif
      ! check if halo component exist
      status = uns_get_range(ident,"halo",n,first,last);
      if (status.gt.0) then ! prepare halo data for saving
         ok=uns_set_array_f(idento,"halo","mass",mass(first)  ,n);
         ok=uns_set_array_f(idento,"halo","pos" ,pos (1,first),n);
         ok=uns_set_array_f(idento,"halo","vel" ,vel (1,first),n);
      endif
      ! check if bulge component exist
      status = uns_get_range(ident,"bulge",n,first,last);
      if (status.gt.0) then ! prepare bulge data for saving
         ok=uns_set_array_f(idento,"bulge","mass",mass(first)  ,n);
         ok=uns_set_array_f(idento,"bulge","pos" ,pos (1,first),n);
         ok=uns_set_array_f(idento,"bulge","vel" ,vel (1,first),n);
      endif

      ! save data
      call uns_save(idento)
      call uns_close_out(idento)

      end








