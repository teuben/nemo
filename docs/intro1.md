# A birds eye view of the stellar dynamics toolbox NEMO

NEMO is a collection of (Unix) programs. Each program is specialized
to do certain things, and you then orchestrate a simulation by scripting a
series of these programs.

## Installing NEMO

An example to use it (*this assumes somebody has installed it for you*) in
the csh shell:

       source /astromake/opt/nemo/git/nemo_start.csh

An example to install it in a bash shell on MacOSX with pgplot
installed via Homebrew:

       curl -O https://teuben.github.io/nemo/install_nemo
       chmod +x install_nemo
       ./install_nemo brew=1
       source nemo/nemo_start.sh

The last **source** command can then be added to your **.cshrc** (csh shell) or
**.bashrc** (linux bash) or **.bash_profile** (mac bash)

## Some properties of NEMO

Your (Unix) shell will now have been modified and a large number
of new commands are available. Much like other Unix commands
these NEMO commands will:
    

      * live in $NEMOBIN, e.g.                    ls $NEMOBIN
      * have a Unix man page, e.g.                man mkplummer
      * have a --help option, e.g.                mkplummer --help
      * have a help= system keyword, e.g.         mkplummer help=\?
        and describe the program keywords, e.g.   mkplummer help=h
      * use Unix pipes                            mkplummer - 10 | tsf -

A few other NEMO things useful to know

      * recompile/add a NEMO program              mknemo ccdsky
      * parameters can have math                  vscale='1/sqrt(2)'
      * lists can have implied do loops           ycol=2:12:2  color=2,3::9,4
      * parameters can refer to others            dtadj='$deltat/2'
      * parameters can be stored in a file        rad=@rad.in
      * nemoinp expression parser                 nemoinp 'c,G,h,pi,iflt(1,2,3,4)'
        1 AU in km                                nemoinp 'p*pi/(3600*180)'
      * The debug= and yapp= system keywords      debug=1   yapp=plot1.ps/vcps
        control debugging level and graphics out  yapp=1/xs yapp=plot1.png/png
      * external potentials dynamically loaded    potname=plummer potpars=0,10,0.5,0.8
      * body transformations dynamically loaded   xvar='sqrt(x*x+y*y)' yvar=vz
      * "snap" "tab", "orb", "ccd" programs       snapplot, tabplot, orbplot, ccdplot
	
## Example 1: The Gentle Collapse of a Plummer (1911) Sphere


Initialize: 1024 stars in a plummer sphere

       mkplummer p1k 1024
       snapscale p1k p1ks vscale="1/sqrt(2)"
       snapplot p1ks
       snapplot p1ks xvar=x yvar=vx

Evolve (pick one, though runbody6 has issues) for about 20 dynamical times

       hackcode1 p1ks run1.out tstop=20 freqout=10 options=mass,phase,phi,acc
       gyrfalcON p1ks run2.out tstop=20 step=0.1    eps=0.05 kmax=7 give=mxvpa
       runbody6  p1ks run3     tcrit=10 deltat=0.1

Plot some properties

       snapplot run1.out nxy=3,3 times=0,1,2,5,10,15,20 yapp=1/xs
       snapmradii run1.out 0.01,0.1:0.9:0.1,0.99 log=t | tabplot - 1 2:12 line=1,1  color=2,3::9,2 yapp=2/xs
       snapvratio run1.out all phi | tabplot - 1 2 ycoord=1 line=1,1

Convert to a "CCD"

       snapgrid run1.out t0.ccd nx=128 ny=128 times=0
       snapgrid run1.out t1.ccd nx=128 ny=128 times=20

Smooth the CCD to 0.1 FWHM 

       ccdsmooth t0.ccd t0s.ccd 0.1
       ccdsmooth t1.ccd t1s.ccd 0.1

Plot the CCD

       ccdplot t0s.ccd contour=0.1,0.3,0.9,1.4 yapp=3/xs
       ccdplot t1s.ccd contour=0.1,0.3,0.9,1.4 yapp=4/xs

Convert to FITS for comparison with observations, or between each other

       ccdfits t0s.ccd t0s.fits
       ccdfits t1s.ccd t1s.fits
       ds9 t0s.fits
       # load t1s.fits in a new frame to compare/blank

Look at outliers via total M/L=1 intensity map (look for Max and Sum*Dx*Dy)

       ccdstat t0s.ccd
       ccdstat t1s.ccd

## Example 2: Cluster Simulations Correctness

When relaxation occurs in small N clusters,
close encounters cannot be followed with softened potentials.
How do we show the correct orbits?

Here we use the data from RAMSES that Chongchong collected as initial conditions
         
       grep  -v e-11 job332out32.in | sed s/^605/603/ |\
	     tabtos - snap1 nbody,time m,pos,vel headline='grep  -v e-11 job332out32.in | sed s/^605/603/'

       snapstat snap1 all=t exact=t
       snapcenter snap1 snap31 m
       snapscale snap31 snap32 1/14332.49 1/0.357757 'sqrt(0.357757/14332.492353)'
       snapstat snap32 all=t exact=t

Example plotting two stars that are interacting  (note this snap2 does not have the correct units)

       hisf snap2
       runbody2 snap2 run4b 605 deltat=0.1 tcrit=100 
       snapplot run4b/OUT3.snap xrange=19:27 yrange=-23:-15 trak=t 'visib=i==287||i==118' 'color=i==287?0.1:0'
       

## Example 3: Exploring a problem with tkrun

Frequently when you break down a particular problem in the analysis,
you have to fine-tune some parameters. On the command line this
can become tedious, and a graphical interface can be handy, even
though you have a manual script. How to put a GUI on your script?
Enter tkrun in NEMO.

       ./do1
       tkrun ./do1

## Other geeky Unix things that NEMO does

   * Uses hierarchy of "Makefile" for installation
   * Uses hierarchy of "Testfile" to do regression/baseline tests
   * Has $NEMO/src (supported) and $NEMO/usr (unsupported)
   * Uses 'mknemo' to recompile a program
   * Can use $NEMO/opt to override bad system libraries

## Useful software that is NEMO friendly

   * ds9
   * glnemo2
   * python  (astropy, APLpy, ....)
   * ImageMagick
   * Montage

## Some other related good stuff

   * ZENO
   * yt
   * AMUSE
   
