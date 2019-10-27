#! /bin/csh -f
#
#   Creates a function, with error bars, which is not quite a straight line
#   then fit a straight line, and get the nominal errors in a and b, 
#   calculates the chi2 map around this best fitted value and shows
#   how by adding in a fake error (q) such that the minimum chi2 is now 1,
#   the +1, +4 and +9 are about the errors in a and b from a classic lsq fit.
#   
#   After this program has run, load the $out.2.fits file into ds9,
#   and overlay the $out.reg file, which is the 1-sigma ellipse 
#   according to fit. Load contourlevels 2,5,10,17, which will be the
#   1,2,3 and 4-sigma chi2 contours.
#  
#   This renormalization is described by Weiner et al, (2001, ApJ 546, 931)
#
#   This script needs NEMO and ds9 for visualization.
#
#   12-nov-2004    --Created--                      Peter Teuben

set n=10
set sig=0.01
set cleanup=1
set apar=0.0:0.5:0.01
set bpar=0.5:1.0:0.01
set out=model1

#                 poor man's command line parser
foreach _a ($*)
  set $_a
end

set tmp=tmp$$

#                 create a table x,y,dy with y(x) some function
nemoinp 1:$n |\
    tabmath - - "%1/$n,sqrt(%2)+rang(0,$sig),$sig" all seed=-1 > $tmp.0

#                 fit a straight line in different ways
tablsqfit $tmp.0 out=$tmp.1 tab=t               > $tmp.fit1
tablsqfit $tmp.0 dycol=3 out=$tmp.2 tab=t       > $tmp.fit2
tabnllsqfit $tmp.0 out=$tmp.3                   > $tmp.fit3
tabnllsqfit $tmp.0 dycol=3 out=$tmp.4           > $tmp.fit4

#                 print out some debugging info
set chi2=`grep ^rms2/ $tmp.fit4 | awk '{print $2}'`
set q=`nemoinp "$sig*sqrt($chi2-1)"`
set sq=`nemoinp "$q*$q+$sig*$sig"`
echo sig=$sig
echo q=$q
echo chi2=$chi2

#                 get the best fit for the line (a=intercept, b=slope)
set a=`grep ^a= $tmp.fit4 | awk '{print $2}'`
set b=`grep ^b= $tmp.fit4 | awk '{print $2}'`
set da=`grep ^a= $tmp.fit4 | awk '{print $3}'`
set db=`grep ^b= $tmp.fit4 | awk '{print $3}'`
echo Linefit: y=a+bx:
echo a=$a  sig_a=$da
echo b=$b  sig_b=$db
#                 create a region file for ds9
echo "linear;ellipse($a,$b,$da,$db,0.0)"      > $out.reg
echo "linear;ellipse($a,$b,0.001,0.001,0.0)" >> $out.reg

#                 for the best fit, show what chi2 and the renormalized chi2+q2 are:
echo -n "sig^2      : " 
tabmath $tmp.0 - "($a+$b*%1-%2)/%3,%4*%4" | tabhist - 5 yapp=/null |& grep Sum | awk -F: '{print $2}'
echo -n "sig^2 + q^2: " 
tabmath $tmp.0 - "($a+$b*%1-%2),%4*%4/$sq" | tabhist - 5 yapp=/null |& grep Sum | awk -F: '{print $2}'

set a=(`nemoinp $apar`)
set b=(`nemoinp $bpar`)
set naxis1=$#a
set naxis2=$#a
set crval1=$a[1]
set crval2=$b[1]
set cdelt1=(`nemoinp $a[2]-$a[1]`)
set cdelt2=(`nemoinp $b[2]-$b[1]`)
set crpix1=1.0
set crpix2=1.0

#                 make a fits header for tabfits, so we get the WCS right
echo "# NAXIS1  =  $naxis1"  > $tmp.fh
echo "# NAXIS2  =  $naxis2" >> $tmp.fh
echo "# CRVAL1  =  $crval1" >> $tmp.fh
echo "# CRVAL2  =  $crval2" >> $tmp.fh
echo "# CDELT1  =  $cdelt1" >> $tmp.fh
echo "# CDELT2  =  $cdelt2" >> $tmp.fh
echo "# CRPIX1  =  $crpix1" >> $tmp.fh
echo "# CRPIX2  =  $crpix2" >> $tmp.fh
echo "# CTYPE1  =  'a'"     >> $tmp.fh
echo "# CTYPE2  =  'b'"     >> $tmp.fh

#                  pre-prend the fits header
cat $tmp.fh > $tmp.5
cat $tmp.fh > $tmp.6

#                  loop over a and b (in the right FITS order) and create the table
foreach b (`nemoinp $bpar`)
 echo -n .
 foreach a (`nemoinp $apar`)
   tabmath $tmp.0 - "($a+$b*%1-%2)/%3,%4*%4" | tabhist - 5 yapp=/null |& grep Sum | awk -F: '{print $2}' >> $tmp.5
   tabmath $tmp.0 - "($a+$b*%1-%2),%4*%4/$sq" | tabhist - 5 yapp=/null |& grep Sum | awk -F: '{print $2}' >> $tmp.6
 end
end

#                    convert the table to a fits file
tabfits $tmp.5  $out.1.fits nx=$naxis1 ny=$naxis2 nz=1
tabfits $tmp.6  $out.2.fits nx=$naxis1 ny=$naxis2 nz=1

echo "Now load image and region:     nds9 $out.2.fits"
echo "ds9->Region->Load Regions...   $out.reg"

if ($cleanup) rm -f $tmp.*
