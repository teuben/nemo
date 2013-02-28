/*
 *
 * CLFIND3 : 
 *    Find clumps in a x-y-v data cube
 *    based on the algorithm described in
 *    Williams, de Geus, & Blitz 1994, ApJ, 428, 693
 *    (please cite this paper in publications that use this procedure)
 *    Re-implemented via IDL from Fortran into C, just to prove a point
 *
 *  Converted from fortran to IDL:          11 Nov 1995  jpw
 *  Complete rewrite using search3d:        29 Mar 2004  jpw
 *  Converted to C in NEMO as CLFIND3:      28 Feb 2013  pjt
 * 
 *
 */

#include <stdinc.h>
#include <strlib.h>
#include <getparam.h>
#include <image.h>
#include <moment.h>
#include <axis.h>
#include <mdarray.h>

string defv[] = {
    "in=???\n             Input file name",
    "out=???\n            Output clump identification file name",
    "step=0.05\n          Contour step",
    "start=1\n            Starting step",
    "VERSION=0.0\n	  28-feb-2013 PJT",
    NULL,
};

string usage = "ClumpFind in 2D or 3D";

string cvsid = "$Id$";

nemo_main()
{
  warning("CLFIND3: testing an IDL to C conversion");
}


#if 0

;...................... START DEFREG ......................
pro defreg,nwork,npix0,reg,nreg
; define regions at contour level nwork
; CLFIND.CB
;-----------------------------------------------------------------------------
common clfindblk1,     infile,levs0,dlevs,ncl,clump_peak
common clfindblk2,     data,assign,nx,ny,nv,bx,by,bv


; get all points at this contour level
levmin=levs0+nwork*dlevs
levmax=levs0+(nwork+1)*dlevs
levpix=where(data ge levmin AND data lt levmax,npix0)
if(npix0 eq 0) then return

; create region array for working level
reg=intarr(nx,ny,nv)
nreg=0

; set all pixels not in level to -1
reg(levpix)=1
reg=reg-1

; loop through finding sucessive unassigned peaks -> new regions
pix1=where(reg eq 0,npix1)
while (npix1 gt 0) do begin
  nreg=nreg+1
  dpeak=max(data(pix1),peak)
  kpeak=pix1(peak)/(nx*ny)
  jpeak=pix1(peak)/nx - kpeak*ny
  ipeak=pix1(peak) - (jpeak+kpeak*ny)*nx
  pix=search3d(data,ipeak,jpeak,kpeak,levmin,levmax,/diagonal)
  reg(pix)=nreg
  pix1=where(reg eq 0,npix1)
endwhile


return
end
;......................  END DEFREG  ......................


;...................... START DEFCLUMP ......................
pro defclump,nwork,reg,nreg,nnew
; define clumps at working contour level
; CLFIND.CB
;-----------------------------------------------------------------------------
common clfindblk1,     infile,levs0,dlevs,ncl,clump_peak
common clfindblk2,     data,assign,nx,ny,nv,bx,by,bv


levmin=levs0+nwork*dlevs
levmax=999999
nnew=0
if (nreg eq 0) then return

; extend each defined region upwards and see
; if it merges with previously defined clumps
for nr=1L,nreg do begin
  pix1=where(reg eq nr,npix1)
  dpeak=max(data(pix1),peak)
  kpeak=pix1(peak)/(nx*ny)
  jpeak=pix1(peak)/nx - kpeak*ny
  ipeak=pix1(peak) - (jpeak+kpeak*ny)*nx
  pix2=search3d(data,ipeak,jpeak,kpeak,levmin,levmax,/diagonal)
  apix2=assign(pix2)

  if(max(apix2) eq 0) then begin
;    print,"Found new clump!"
    nnew=nnew+1
    assign(pix2)=ncl+nnew
    dpeak=max(data(pix2),peak)
    kpeak=pix2(peak)/(nx*ny)
    jpeak=pix2(peak)/nx - kpeak*ny
    ipeak=pix2(peak) - (jpeak+kpeak*ny)*nx
    clump_peak(ncl+nnew,0)=ipeak
    clump_peak(ncl+nnew,1)=jpeak
    clump_peak(ncl+nnew,2)=kpeak

  endif else begin
;   define list of merged clumps
    nc=apix2(uniq(apix2,sort(apix2)))
    if(min(nc) eq 0) then begin

      if(n_elements(nc) eq 2) then begin
;       just one clump above working level, so extend it
;       print,"Extending clump ",nc(1)
        assign(pix2)=nc(1)

      endif else begin
;       more than one clump so this region is a merger
        clump_merge=nc(1:*)
;        print,"Merging clumps ",clump_merge
        ic=clump_peak(clump_merge,0)
        jc=clump_peak(clump_merge,1)
        kc=clump_peak(clump_merge,2)

;       go through pixel by pixel and assign to nearest clump
        for nr1=0L,npix1-1 do begin
          k=pix1(nr1)/(nx*ny)
          j=pix1(nr1)/nx - k*ny
          i=pix1(nr1) - (j+k*ny)*nx
          d=(i-ic)^2 + (j-jc)^2 + (k-kc)^2
          dmin=min(d,m)
          assign(pix1(nr1))=clump_merge(m)
        endfor
      endelse
    endif

  endelse
endfor


return
end
;......................  END DEFCLUMP  ......................


;...................... START TESTBAD ......................
pro testbad,nmin,nbad
; sorts clumps in order of peak flux and
; rejects those with number of pixels <= nmin
; CLFIND.CB
;-----------------------------------------------------------------------------
common clfindblk1,     infile,levs0,dlevs,ncl,clump_peak
common clfindblk2,     data,assign,nx,ny,nv,bx,by,bv

nc=indgen(ncl)+1
dmax=data(clump_peak(nc,0),clump_peak(nc,1),clump_peak(nc,2))
new_order=reverse(sort(dmax))

ncl_new=0
nbad=0
assign0=assign

for n1=1L,ncl do begin
  n0=new_order(n1-1)+1
  iclp=where(assign0 eq n0,count)

  if (count le nmin) then begin
    nbad=nbad+1
    assign(iclp)=0
  endif else begin
    ncl_new=ncl_new+1
    assign(iclp)=ncl_new
  endelse

endfor


return
end
;......................  END TESTBAD  ......................



;...................... START MK_HDR ......................
pro mk_hdr,header,header_out
;------------------------------------
; create fits header for output file
; (basically a copy of input header
; with some program information)
;------------------------------------
; CLFIND.CB
;-----------------------------------------------------------------------------
common clfindblk1,     infile,levs0,dlevs,ncl,clump_peak
common clfindblk2,     data,assign,nx,ny,nv,bx,by,bv


s=size(header)
nlines=s(1)
header_out=strarr(25)
header_out(0)="SIMPLE  =                    T /"
header_out(1)="BITPIX  =                   16 /"
header_out(2)="NAXIS   =                    3 /"
n=2
for i=0L,nlines-1 do begin
  if (strpos(header(i),'NAXIS1') eq 0) then begin
    n=n+1 & header_out(n)=header(i)
  endif
  if (strpos(header(i),'NAXIS2') eq 0) then begin
    n=n+1 & header_out(n)=header(i)
  endif
  if (strpos(header(i),'NAXIS3') eq 0) then begin
    n=n+1 & header_out(n)=header(i)
  endif
  if (strpos(header(i),'CRVAL1') eq 0) then begin
    n=n+1 & header_out(n)=header(i)
  endif
  if (strpos(header(i),'CRVAL2') eq 0) then begin
    n=n+1 & header_out(n)=header(i)
  endif
  if (strpos(header(i),'CRVAL3') eq 0) then begin
    n=n+1 & header_out(n)=header(i)
  endif
  if (strpos(header(i),'CRPIX1') eq 0) then begin
    n=n+1 & header_out(n)=header(i)
  endif
  if (strpos(header(i),'CRPIX2') eq 0) then begin
    n=n+1 & header_out(n)=header(i)
  endif
  if (strpos(header(i),'CRPIX3') eq 0) then begin
    n=n+1 & header_out(n)=header(i)
  endif
  if (strpos(header(i),'CDELT1') eq 0) then begin
    n=n+1 & header_out(n)=header(i)
  endif
  if (strpos(header(i),'CDELT2') eq 0) then begin
    n=n+1 & header_out(n)=header(i)
  endif
  if (strpos(header(i),'CDELT3') eq 0) then begin
    n=n+1 & header_out(n)=header(i)
  endif
  if (strpos(header(i),'CTYPE1') eq 0) then begin
    n=n+1 & header_out(n)=header(i)
  endif
  if (strpos(header(i),'CTYPE2') eq 0) then begin
    n=n+1 & header_out(n)=header(i)
  endif
  if (strpos(header(i),'CTYPE3') eq 0) then begin
    n=n+1 & header_out(n)=header(i)
  endif
end

header_out(n+1)=string(levs0,format='("CLFIND: lowest contour level  =",f6.2," /")')
header_out(n+2)=string(dlevs,format='("CLFIND: contour increment =",f6.2," /")')
header_out(n+3)="END"+string(replicate(32b,77))

; get rid of blank lines (if any) at end
header_out=header_out(0:n+3)

return
end
;......................  END MK_HDR  ......................

;...................... START CLFIND ......................
pro clfind,file=file,low=low,inc=inc,log=log

forward_function defreg, defclump, testbad, mk_hdr

; read in fits cube and call clump finding routines

; CLFIND.CB
;-----------------------------------------------------------------------------
common clfindblk1,     infile,levs0,dlevs,ncl,clump_peak
common clfindblk2,     data,assign,nx,ny,nv,bx,by,bv

;-----------------------------------------------------------------------------
common header_block,    naxis1,crpix1,cdelt1,crval1,ctype1,cd1_1, $
                        naxis2,crpix2,cdelt2,crval2,ctype2,cd2_2, $
                        naxis3,crpix3,cdelt3,crval3,ctype3
;-----------------------------------------------------------------------------

if NOT keyword_set(file) then begin
  print,'PROCEDURE clfind,file=filename,low=low,inc=inc,[/log]'
  print,'------------------------------------------------------------'
  print,'filename = name of the FITS data cube (in quotes)'
  print,'low      = lowest contour level'
  print,'inc      = contour level increment'
  print,'/log       for screen output copied to clfind.log'
  print,'------------------------------------------------------------'
  return
endif 

print,"----------------------------------------------------------------"
print,"CLFIND: ",systime()
print,"----------------------------------------------------------------"
print,format='("Filename = ",a0)',file
print,format='("Lowest contour level =",f6.2)',low
print,format='("Contour increment =",f6.2)',inc
print,"----------------------------------------------------------------"
if keyword_set(log) then begin
  openw,1,'clfind.log'
  printf,1,"----------------------------------------------------------------"
  printf,1,"CLFIND: ",systime()
  printf,1,"----------------------------------------------------------------"
  printf,1,format='("Filename = ",a0)',file
  printf,1,format='("Lowest contour level =",f6.2)',low
  printf,1,format='("Contour increment =",f6.2)',inc
  printf,1,"----------------------------------------------------------------"
endif

; read in data + header
data=readfits(file+'.fits',header,/silent)
readhd,header
nx=long(naxis1)
ny=long(naxis2)
nv=long(naxis3)

; save for common block
infile=file
levs0=low
dlevs=inc

; mask out bad data (NaN replaced by -999.9)
bad=where(finite(data) eq 0,count)
if (count gt 0) then data(bad)=-999.9

; initialize clump assignment arrays
max_clumps=9999
assign=intarr(nx,ny,nv)
clump_peak=intarr(max_clumps,3)


; MAIN CLUMP FINDING LOOP
ncl=0
t0=systime(1)
nlev=1+fix((max(data)-levs0)/dlevs)

for nwork=nlev-1,0,-1 do begin
  defreg,nwork,npix,reg,nreg
  print,format='("Contour level ",f6.2,": ",i5," pixels ",i5," regions ",$)'$
       ,levs0+nwork*dlevs,npix,nreg
  if keyword_set(log) then $
    printf,1,format='("Contour level ",f6.2,": ",i5," pixels ",i5," regions ",$)'$
            ,levs0+nwork*dlevs,npix,nreg

  defclump,nwork,reg,nreg,nnew
  print,format='(i5," new clumps")',nnew
  if keyword_set(log) then printf,1,format='(i5," new clumps")',nnew
  ncl=ncl+nnew

endfor

; reject clumps with fewer than npixmin pixels
npixmin=5
testbad,npixmin,nbad
ncl=ncl-nbad
print,format='(i4," clumps found (",i3," rejected)")',ncl,nbad
print,"================================================================"
if keyword_set(log) then begin
  printf,1,format='(i4," clumps found (",i3," rejected)")',ncl,nbad
  printf,1,"================================================================"
endif

; write assign file to a fits file
outfile=infile+'.fits.clf'
mk_hdr,header,outheader
writefits,outfile,assign,outheader
print,"Writing output file: ",outfile
print,""
print,"CLFIND exits sucessfully"

delta_t=(systime(1)-t0)/60.0
print,format='(f5.1," minutes elapsed")',delta_t

if keyword_set(log) then begin
  printf,1,format='("Writing output file: ",a0)',outfile
  printf,1,format='(f5.1," minutes elapsed")',delta_t
  close,1
endif

return
end
;......................  END CLFIND  ......................

#endif
