include $(NEMOLIB)/makedefs
#
#FC = g77
#FC = pgf77
#
#FFLAGS = -O3
#FFLAGS = -ffortran-bounds-check -O3
#FFLAGS = -O4 -tp p7 -pc 64
#

SOURCE = nbody4.f \
bindat.f binout.f binpop.f block.f bodies.f centre.f \
clint.f cloud.f cloud0.f cmbody.f cmcorr.f coal.f comenv.f core.f \
corerd.f cputim.f data.f define.f degen.f delay.f derqp3.f \
derqp4.f dfcm2.f dfsp.f dfsp2.f dgcore.f difsy3.f difsy4.f \
dtchck.f efac2.f efac3.f endreg.f erel3.f erel4.f escape.f events.f \
evolve.f expand.f expel.f extend.f fcloud.f fcorr.f fdisk.f ffdot.f \
ffdot2.f ffdot3.f fhalo.f findj.f flush.f flyby.f fnuc.f fpert.f \
fpoly1.f freeze.f gcinit.f gcint.f giant.f giant3.f global.f \
gntage.f hcorr.f hiarch.f hicirc.f hidat.f himax.f himod.f hipop.f \
hirect.f histab.f hivel.f hmdot.f hmdot2.f hotsys.f \
hrdiag.f hrplot.f iblock.f ichain.f imf.f imfbd.f impact.f impuls.f \
input.f jacobi.f insert.f instar.f intev.f intide.f kepler.f kick.f \
ksapo.f kscorr.f ksin2.f ksinit.f ksint.f kslist.f ksmod.f ksperi.f \
kspert.f kspoly.f kspred.f ksrect.f ksreg.f ksres.f ksres2.f ksterm.f \
kstide.f kstran.f lagr.f lagr2.f levels.f matrix.f mdot.f mdot2.f \
merge.f merge2.f mix.f mloss.f mlwind.f modify.f mrenv.f mtrace.f \
mturn.f mydump.f nbpot.f nbrem.f nbrest.f nclist.f newreg.f newsys.f \
nlmod.f ntint.f offset.f orbit.f output.f peri.f period.f \
permit.f pfac.f phidot.f potn0.f qpmod3.f qpmod4.f quad.f ran2.f \
rchain.f reflct.f remove.f rename.f reset.f reset2.f resolv.f \
rkint.f rl.f rpmax.f rsort.f scale.f setup.f shock.f short.f \
shrink.f sort1.f stability.f stabl3.f stabl4.f stablz.f \
star.f start3.f start4.f status.f stumpf.f subsys.f tail0.f \
tcirc.f tides.f tides2.f tnew.f touch.f tperi.f tpert.f trans3.f \
trans4.f trdot.f trflow.f triple.f tstab.f tstep.f units.f \
unpert.f update.f verify.f xtrnl0.f xtrnld.f xtrnlf.f xtrnlp.f \
xtrnlt.f xtrnlv.f xvpred.f zare.f zcnsts.f zero.f zfuncs.f \
adjust.f check.f energy.f fpoly0.f fpolyi.f g6acc.f \
inext.f intgrt.f nbcorr.f nblist.f pipes.f poti.f potn2.f \
search.f start.f stepi.f stepk.f steps.f subint.f \
cfuncs.f chain.f chstab.f const.f cstab2.f cstab3.f cstab4.f cstab5.f \
derqp.f difsy1.f erel.f hpsort.f inclin.f invert.f ksphys.f physks.f \
qforce.f qpmod.f r2sort.f recoil.f redraw.f select.f slow.f stablc.f \
swcond.f switch.f touch2.f transk.f transq.f transx.f vector.f xtf.f \
xtrnlu.f ycopy.f ysave.f \
absorb.f chdata.f chf.f chfind.f chinit.f chlist.f chmod.f chpert.f \
chpot.f chterm.f cmlist.f expel2.f fchain.f ghost.f kcpert.f reduce.f \
reinit.f renew.f setsys.f tchain.f xcpred.f xtpert.f gpfree.f gpsend.f

OBJECTS = $(SOURCE:.f=.o)

nbody4:	$(OBJECTS)
	$(FC) $(FFLAGS) $(OBJECTS) -o nbody4

clean:
	\rm -f *.o nbody4

install: nbody4
	cp $? $(NEMOBIN)
