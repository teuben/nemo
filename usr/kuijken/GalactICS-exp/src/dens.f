
      function dens(r,z)
      parameter (pi=3.1415926535)
      common /flags/ idiskflag, ibulgeflag, ihaloflag
     
      if( idiskflag .eq. 1 ) then
         addens = appdiskdens(r,z)
         dens=totdens(r,z)
         dens = dens - addens
      else
         dens=totdens(r,z)
      endif
      return
      end
