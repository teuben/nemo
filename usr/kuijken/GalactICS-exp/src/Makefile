#Makefile for particle generators of the combined df - June 8/94 - JJD
#
.f.o:
	$(F77) $(FFLAGS) -c $*.f

F77=f77
CC=cc
CFLAGS = -O -DRINGASCII -DASCII
LIBS = -lm
#These flags for SunOS:
FFLAGS =  -e -C -O 
FLIBS = -L/usr/X11/lib/ -lX11 -lpgplot 

#Replace FFLAGS and FLIBS above with the following for linux:
#FFLAGS = -ffixed-line-length-0 -O
#FLIBS  = -L/usr/X11R6/lib/ -lX11 -lpgplot 
#

ALLFILES= dbh genhalo genbulge gendisk diskdf tobinary \
     getfreqs dbh vcirc plothalo plotforce 

#FFLAGS = -O -extend_source
#CFLAGS = -O2 -DRINGASCII -DASCII -DTP
#LIBS = -lm -lots -lfor


all: $(ALLFILES)

dbh: dbh.o dbhplot.o polardens.o bulgepotential.o totdens.o halopotential.o pot.o diskdens.o dens.o appdiskpot.o plgndr1.o bulgedenspsi.o densrpsi.o polarbulgedens.o polarhalodens.o appdiskdens.o halodens.o erfcen.o modstamp.o
	$(F77) $(FFLAGS) $(FLIBS) dbh.o dbhplot.o polardens.o bulgepotential.o totdens.o halopotential.o pot.o diskdens.o dens.o appdiskpot.o plgndr1.o bulgedenspsi.o densrpsi.o polarbulgedens.o polarhalodens.o appdiskdens.o halodens.o erfcen.o modstamp.o  -o dbh

genhalo: genhalo.o readmassrad.o query.o ran1.o readharmfile.o pot.o halodens.o appdiskpot.o plgndr1.o densrpsi.o erfcen.o
	$(F77) $(CFLAGS) genhalo.o readmassrad.o query.o ran1.o readharmfile.o appdiskpot.o plgndr1.o densrpsi.o pot.o halodens.o erfcen.o -o genhalo

genbulge: genbulge.o readmassrad.o query.o ran1.o readharmfile.o pot.o bulgedens.o appdiskpot.o plgndr1.o bulgedenspsi.o erfcen.o
	$(F77) $(CFLAGS) genbulge.o readmassrad.o query.o ran1.o readharmfile.o appdiskpot.o plgndr1.o bulgedenspsi.o pot.o bulgedens.o erfcen.o -o genbulge

gendisk: gendisk.o ran1.o query.o \
	readdiskdf.o golden.o simpson.o diskdf5ez.o diskdensf.o readharmfile.o sigr2.o sigz2.o omekap.o diskdens.o splined.o splintd.o diskdf3ez.o diskdenspsi.o pot.o fnamidden.o appdiskpot.o plgndr1.o rcirc.o
	$(F77) $(CFLAGS) gendisk.o ran1.o query.o diskdf5ez.o diskdensf.o readharmfile.o sigr2.o sigz2.o omekap.o diskdens.o splined.o splintd.o diskdf3ez.o diskdenspsi.o pot.o fnamidden.o appdiskpot.o plgndr1.o rcirc.o \
	readdiskdf.o golden.o simpson.o -o gendisk

diskdf: diskdf.o diskdf5intez.o diskdensf.o splintd.o splined.o readharmfile.o sigr2.o sigz2.o omekap.o pot.o diskdf3intez.o appdiskpot.o plgndr1.o fnamidden.o rcirc.o diskdens.o modstamp.o
	$(F77) $(FFLAGS)  $(FLIBS) diskdf.o diskdf5intez.o diskdf3intez.o appdiskpot.o plgndr1.o diskdensf.o splintd.o splined.o readharmfile.o sigr2.o sigz2.o fnamidden.o rcirc.o diskdens.o omekap.o pot.o modstamp.o  -o diskdf

getfreqs: getfreqs.o pot.o appdiskpot.o plgndr1.o erfcen.o
	$(F77) $(FFLAGS)  $(FLIBS) getfreqs.o pot.o appdiskpot.o plgndr1.o erfcen.o -o getfreqs

tobinary: tobinary.o
	$(F77) $(CFLAGS) tobinary.o -o tobinary

vcirc: vcirc.o query.o readharmfile.o force.o appdiskforce.o appdiskpot.o plgndr1.o
	$(F77) $(CFLAGS) vcirc.o query.o  readharmfile.o force.o appdiskforce.o appdiskpot.o plgndr1.o -o vcirc

plothalo: plothalo.o contourden.o readdiskdf.o diskdens.o halodens.o bulgedens.o bulgedenspsi.o diskdensf.o readharmfile.o sigr2.o densrpsi.o omekap.o pot.o splined.o appdiskpot.o splintd.o diskdenspsi.o erfcen.o plgndr1.o  modstamp.o
	$(F77) $(FFLAGS) $(FLIBS)  plothalo.o contourden.o readdiskdf.o diskdens.o halodens.o diskdensf.o bulgedens.o bulgedenspsi.o readharmfile.o sigr2.o densrpsi.o omekap.o pot.o splined.o appdiskpot.o splintd.o diskdenspsi.o erfcen.o plgndr1.o modstamp.o  -o plothalo

plotforce: plotforce.o readharmfile.o readdiskdf.o force.o appdiskforce.o appdiskpot.o diskdensf.o sigr2.o omekap.o plgndr1.o diskdens.o splined.o splintd.o diskdenspsi.o pot.o modstamp.o
	$(F77) $(FFLAGS) $(FLIBS)  plotforce.o readharmfile.o readdiskdf.o force.o appdiskforce.o appdiskpot.o diskdensf.o sigr2.o omekap.o plgndr1.o diskdens.o splined.o splintd.o diskdenspsi.o pot.o  modstamp.o  -o plotforce

#clean:
#       rm -f *.o genhalo genbulge gendisk genglob genrings \
#         diskdf getfreqs idz mring dbh vcirc tobinary spherden orbit mutpot
#       (cd ../bin; rm -f *.o genhalo genbulge gendisk genglob genrings \
#         diskdf getfreqs idz mring dbh vcirc tobinary spherden)
#install:
#       cp dbh getfreqs diskdf gendisk genbulge genhalo genglob \
#       vcirc tobinary spherden mergerv ../bin


clean:
	rm -f *.o $(ALLFILES)
	(cd ../bin; rm -f *.o $(ALLFILES))

install:
	cp -f $(ALLFILES) mergerv ../bin


readmassrad.o: readmassrad.c
simpson.o: simpson.c
ran1.o: ran1.c
query.o: query.c
golden.o: golden.c
tobinary.o: tobinary.c

