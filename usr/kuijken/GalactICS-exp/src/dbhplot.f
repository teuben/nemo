      subroutine dbhplot(apot,lmax,nr,dr)
      parameter(ipmax=1000)
      parameter(jpmax=20000)
      integer nrp,iplot(ipmax)
      real apot(20,0:jpmax),rp(ipmax),pplot(ipmax)
      data ifirst /0/
      save rp,ifirst,nrp
      character*40 toplbl

      write(*,*) 'DBHPLOT: ',lmax,nr,dr,ifirst
      if (ifirst.eq.0) then
         do ir=1,ipmax
            iplot(ir)=0
         enddo
         ifirst=1
         nrp=1
         rp(nrp)=log10(dr)
         iplot(nrp)=1 
         write(*,*) 'DBHPLOT=> ',1,rp(nrp),nrp
         do ir=1,nr
            r=log10(ir*dr)
            if (r-rp(nrp).gt.0.01) then
               nrp=nrp+1
               rp(nrp)=r
               iplot(nrp)=ir
               write(*,*) 'DBHPLOT=> ',ir,r,nrp
               if (nrp.ge.ipmax) then
                  write(*,*) 'Not able to plot all points in dbhplot.'
                  write(*,*) 'Increase ipmax parameter & recompile!'
                  goto 4
               endif
               if (ir.ge.jpmax) then
                  write(*,*) 'ir=',ir,' too big for ',jpmax
                  stop
               endif
            endif
         enddo
         write(*,*) '  dbhplot::nrp=',nrp
      endif

c     big cheat, since in the code below on 64-bit computers there
c     is some array indexing wrong
      return

 4    do l=0,lmax,2
         j=l/2+1
         if (l.eq.0) then
            pmax=0.
         else
            pmax=-1e31
         endif
         pmin=1e31
         do irp=1,nrp
            pp=apot(j,iplot(irp))    
            if (l.eq.0) then
               pp=min(pp,0.)
            else
               pmax=max(pmax,pp)
            endif
            pmin=min(pmin,pp)
            pplot(irp)=pp
         enddo
         call pgenv(log10(dr),log10(nr*dr),pmin,pmax,0,1)
         call pgline(nrp,rp,pplot)
         write(toplbl,'(''l='',i3)') l
         call pglabel('log10(R)','POT',toplbl)
      enddo
      return
      end
