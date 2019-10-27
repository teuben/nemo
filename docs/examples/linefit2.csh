#! /bin/csh -f
#
#
#  performance:
#  - 1000 calls to this script with plotting takes 3.5 mins (cpu 50%)
#  - no plotting takes 12"
#  - removing lablsqfit as well takes 10"

set x=7:11.5:0.05      # the X (log(X) really) coordinates
set sx=0.5             # error in X
set sy=0.5             # error in Y
set a=4                # intercept
set b=0.55             # slope
set seed=-1            # seed
set doplot=1           # plotting?

foreach arg ($*)
  set $arg
end

set tab=tmp.linefit2.tab

set xr=(xmin=6 xmax=12)
set yr=(ymin=6 ymax=12)
set yr2=(ymin=-6 ymax=6)


# col 1 = x (no error)
#     2 = y (no error)
#     3 = x' = x+e_x (observed)
#     4 = y' = y+e_y (observed)
#     5 = y'/x'            (the wrong way)
#     6 = log(10^y'/10^x') (the correct way)
nemoinp $x | tabmath - - "$a+$b*%1,%1+rang(0,$sx),%2+rang(0,$sy),%4/%3,log(10**%4/10**%3)" seed=$seed > $tab

if ($doplot) then
  tabplot $tab 3 4 point=2,0.2 yapp=1/xs   $xr $yr
  tabplot $tab 1 4 point=2,0.2 yapp=2/xs   $xr $yr
  tabplot $tab 3 6 point=2,0.2 yapp=3/xs   $xr $yr2
endif

tablsqfit $tab 3 4
tablsqfit $tab 1 4
tablsqfit $tab 3 6

linreg $tab 3 4
linreg $tab 1 4
linreg $tab 3 6
