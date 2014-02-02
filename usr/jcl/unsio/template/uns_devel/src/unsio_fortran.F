! =====================================================================
! Copyright Jean-Charles LAMBERT - 2008-2014
!           Centre de donneeS Astrophysiques de Marseille (CeSAM)
! e-mail:   Jean-Charles.Lambert@lam.fr                                      
! address:  Aix Marseille Universite, CNRS, LAM        
!           Laboratoire d'Astrophysique de Marseille                          
!           Pole de l'Etoile, site de Chateau-Gombert                         
!           38, rue Frederic Joliot-Curie                                     
!           13388 Marseille cedex 13 France                                   
!           CNRS UMR 7326                                                     
! =====================================================================
! The following program shows how to use UNSIO library from
! Fortran program
!
! This program reads an unsio compatible snapshot from the command line
! and save it in gadget2 format
!
! Syntaxe : unsio_fortran myinput myoutput select_comp select_time
!
! myinput     -> an unsio compatible input snapshot
! myoutput    -> output filename
! select_comp -> the component to be saved
! select_time -> a range of time which select input snapshot
!
! For more information about how to use UNSIO, visit:
! http://projets.oamp.fr/projects/unsio/
!
! =====================================================================
      program unsio_fortran
      implicit none

      integer lnblnk, narg, uns_init, uns_load, valid, ident,
     $     nbody, cpt, uns_get_value_i,
     $     status

      character ain * 80, aout * 80, atyp * 80, acomp * 80, atim * 80
      character out * 80, ii * 6, simdir * 100

      ! check command line
      narg = command_argument_count() ! #parameters
      if (narg.ne.4) then
          call getarg(0,ain)
         write(0,*) "You must give 4 parameters:"
         write(0,*) ain(1:lnblnk(ain)),
     $        " inputfile outputfile selected_component",
     $        " selected_time"
          write(0,*) ""
          write(0,*) "Example : "

          write(0,*) ain(1:lnblnk(ain)),
     $         " sgs029 sgs029.gas gas  all"
         stop
      endif

      ! read arguments from the command line
      call get_command_argument(1,ain)     ! input snapshot
       
      call get_command_argument(2,aout)    ! output filename

      call get_command_argument(3,acomp)   ! selected component

      call get_command_argument(4,atim)    ! selected time

      ! output file format is gadget2
      atyp = "gadget2"  

      ! ************************
      ! initialyze UNSIO engine
      ! ************************
      ident=uns_init(ain, acomp, atim) ! return identifier for the snaphot

      if (ident.gt.0) then  ! ident must be positive
         valid = 1
         cpt = 0
         call uns_sim_dir(ident,simdir)
         write (0,*) "Simulation DIR = [",simdir(1:lnblnk(simdir)),"]"

         do while (valid .gt. 0) ! loop on all time steps
            valid = uns_load(ident) ! load data belonging to ident snapshot
            write(0,*) "uns_load return=>",valid
            if (valid .gt. 0) then  ! it goes fine 

               !! build output file name
               write (ii,'(I6)') cpt
               out = trim(aout)//"."//trim(adjustl(ii))
               write (0,*) out

               ! get nbody, used later for automatic memory 
               ! allocation in start subroutine
               status=uns_get_value_i(ident,"nsel", nbody) ! get #bodies
               write (0,*) "nbody =", nbody

               !! start main subroutine
               call start(ident,nbody,out,atyp)
               cpt = cpt+1
            endif
         enddo
      endif
      end
      
      ! =====================================================================
      ! START subroutine, read snapshot
      ! =====================================================================
      subroutine start(ident,nbody,out,atyp)
      implicit none
!     input parameters
      integer ident,nbody
!     UNS variable
      integer status
      integer  uns_get_range, uns_save_init
      integer  uns_get_value_f, uns_get_array_f, nsel
      integer  uns_set_array_f, uns_set_value_f, ok
      real *4 time, pos(3,nbody), vel(3,nbody), mass(nbody), rho(nbody),
     $     hsml(nbody), u(nbody)
      character out * 80, atyp * 80, fstructure * 90, interface * 90,
     $     fname * 90
!     various      
      integer n,first,last, idento

! initialyze uns output
      idento=uns_save_init(out,atyp)

      status=uns_get_value_f(ident,"time",time      )
      write (0,*) 'time =',time

      
! prepare time for saving
      ok=uns_set_value_f(idento,"time", time)

! get pos vel mass for all particles

      nsel=uns_get_array_f(ident,"all","pos" ,pos ,nbody)  ! pos
      nsel=uns_get_array_f(ident,"all","vel" ,vel ,nbody)  ! vel
      nsel=uns_get_array_f(ident,"all","mass",mass,nbody)  ! mass
      write (0,*) 'nsel =',nsel,' nbody=',nbody

      call uns_get_file_structure(ident,fstructure)
      write (0,*) "File structrure = [",
     $     fstructure(1:lnblnk(fstructure)),"]"

      call uns_get_interface_type(ident,interface)
      write (0,*) "Interface type  = [",
     $     interface(1:lnblnk(interface)),"]"

      call uns_get_file_name(ident,fname)
      write (0,*) "File name       = [",fname(1:lnblnk(fname)),"]"

      if (fstructure(1:lnblnk(fstructure)).eq."component" ) then
         write (0,*) "it's a component snapshot !!!!"
         
! check if gas component exist
         status = uns_get_range(ident,"gas",n,first,last);
         if (status.gt.0) then  ! prepare gas data for saving
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
         if (status.gt.0) then  ! prepare disk data for saving
            ok=uns_set_array_f(idento,"disk","mass",mass(first)  ,n);
            ok=uns_set_array_f(idento,"disk","pos" ,pos (1,first),n);
            ok=uns_set_array_f(idento,"disk","vel" ,vel (1,first),n);
         endif
! check if stars component exist
         status = uns_get_range(ident,"stars",n,first,last);
         if (status.gt.0) then  ! prepare stars data for saving
            ok=uns_set_array_f(idento,"stars","mass",mass(first)  ,n);
            ok=uns_set_array_f(idento,"stars","pos" ,pos (1,first),n);
            ok=uns_set_array_f(idento,"stars","vel" ,vel (1,first),n);
         endif
! check if halo component exist
         status = uns_get_range(ident,"halo",n,first,last);
         if (status.gt.0) then  ! prepare halo data for saving
            ok=uns_set_array_f(idento,"halo","mass",mass(first)  ,n);
            ok=uns_set_array_f(idento,"halo","pos" ,pos (1,first),n);
            ok=uns_set_array_f(idento,"halo","vel" ,vel (1,first),n);
         endif
! check if bulge component exist
         status = uns_get_range(ident,"bulge",n,first,last);
         if (status.gt.0) then  ! prepare bulge data for saving
            ok=uns_set_array_f(idento,"bulge","mass",mass(first)  ,n);
            ok=uns_set_array_f(idento,"bulge","pos" ,pos (1,first),n);
            ok=uns_set_array_f(idento,"bulge","vel" ,vel (1,first),n);
         endif

      else
         ! it's RANGE type snapshot
         write (0,*) "it's a RANGE type snapshot !!!!"      
         ! we are going to save all data (pos,vel,mass) in the HALO field
         ok=uns_set_array_f(idento,"halo","mass",mass,nsel);
         ok=uns_set_array_f(idento,"halo","pos" ,pos ,nsel);
         ok=uns_set_array_f(idento,"halo","vel" ,vel ,nsel);
      endif
! save data
      call uns_save(idento)

      end








