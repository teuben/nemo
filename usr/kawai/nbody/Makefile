# intel:   icpc/icc
CXX = g++
CC  = gcc
#CXX = icpc
#CC  = icc

# PGINCLUDE =  -I$(PGPLOT_DIR)
G5INCLUDE = -I/home/atsushi/pub/include

#
# Add -DNOGRAPHICS for installation without PGPLOT
#
# following option might increase the performance:
# -fp_reorder, -misalign, -nocleanup, -noinline, -tune,
# -assume trusted_short_alignment, -unroll
#
# definition to generate some variants
# -pthread -DPAR_TRAV: performs tree traverse and force calculation in parallel
# -DP2M2: use Pseudo-particle Multipole Method
#
# direct summation
# CXXFLAGS = -DNOGRAPHICS $(PGINCLUDE) -O4 -tune ev5 -arch ev5
#
# treecode
# CXXFLAGS = -DNOGRAPHICS -DTREE $(PGINCLUDE) -O4 -tune ev5 -arch ev5
# -pthread -DPAR_TRAV
#
# P2M2 treecode
# CXXFLAGS = -DNOGRAPHICS -DTREE $(PGINCLUDE) -O4 -tune ev5 -arch ev5 -DP2M2
# -DUSE_P2M2_ALWAYS
#
# only for test
# CXXFLAGS = -DNOGRAPHICS -DTREE $(PGINCLUDE) -O4 -tune ev5 -arch ev5 -DP2M2 -DTESTCODE -DUSE_P2M2_ALWAYS
#
# for g++ & gcc
# CXXFLAGS = -DNOGRAPHICS $(PGINCLUDE) -O3 -ffast-math -funroll-loops
# CXXFLAGS = -DNOGRAPHICS -DTREE $(PGINCLUDE) -O3 -ffast-math -funroll-loops -DPAR_TRAV -pthread
CXXFLAGS = -DNOGRAPHICS -DTREE $(PGINCLUDE) -O3 -DP2M2
#CXXFLAGS += -ffast-math -funroll-loops 
#CXXFLAGS += -ipo -xW -tpp7 -prefetch

# IONEMO = -L$(NEMOLIB) -lio_nemo -lnemo
# PGLIB     =  -L$(PGPLOT_DIR) -O4 -lcpgplot -lpgplot -lX11 -lUfor -lfor -lFutil -lm -lprof1 -lpdf

G5LIB =  -L/home/atsushi/test1.5/g5pkg/lib -lg5a -lphibdma
G4LIB =  -L/home/atsushi/test1.5/g5pkg/lib -lg4a -lrt
H3LIB = -L/home/atsushi/test1.5/g5pkg/lib -lharp3  -lrt

LIBOBJS = pgetopt.o gravity.o BHtree.o second.o p2m2.o
LIBOBJSG4 = pgetopt.o  gravity_g4.o BHtree_g4.o second.o p2m2_g4.o
LIBOBJSG5 = pgetopt.o  gravity_g5.o BHtree_g5.o second.o p2m2_g5.o

EXPORT_FILES = BHtree.h nbody_particle.h vector.h \
               BHtree.C	gravity.C nbody.C pgetopt.C \
	       second.c design.c p2m2.C p2m2.h library.C\
	       Makefile Makefile.template nbody_guide.tex\
               README COPYRIGHT INSTALL
#               samplein sampleparm\
#	 samplelog samplelog_g4 samplelog_g5

gravity.o : gravity.C  nbody_particle.h
	$(CXX)  $(CXXFLAGS) -c gravity.C
gravity_g4.o : gravity.C  nbody_particle.h
	$(CXX)  -DHARP3 $(CXXFLAGS)  -c gravity.C
	mv gravity.o gravity_g4.o
gravity_g5.o : gravity.C  nbody_particle.h
	$(CXX)  -DGRAPE5 $(G5INCLUDE) $(CXXFLAGS)  -c gravity.C
	mv gravity.o gravity_g5.o

p2m2.o : p2m2.C
	$(CXX)  $(CXXFLAGS) -c p2m2.C

p2m2_g4.o : p2m2.C
	$(CXX)  -DHARP3 $(CXXFLAGS) -c p2m2.C
	mv p2m2.o p2m2_g4.o

p2m2_g5.o : p2m2.C
	$(CXX)  -DGRAPE5 $(G5INCLUDE) $(CXXFLAGS) -c p2m2.C
	mv p2m2.o p2m2_g5.o

libnbody.a:	library.o $(LIBOBJS)
	ar ruv libnbody.a library.o $(LIBOBJS)
	ranlib libnbody.a

libnbodyg5.a:	library_g5.o $(LIBOBJSG5)
	ar ruv libnbodyg5.a library_g5.o $(LIBOBJSG5)
	ranlib libnbodyg5.a

libnbodyg4.a:	library_g4.o $(LIBOBJSG4)
	ar ruv libnbodyg4.a library_g4.o $(LIBOBJSG4)
	ranlib libnbodyg4.a

library.o:	library.C nbody_particle.h
	$(CXX) $(CXXFLAGS) -o library.o -c library.C

library_g5.o:	library.C nbody_particle.h
	$(CXX) -DGRAPE5 $(G5INCLUDE) $(CXXFLAGS) -o library_g5.o -c library.C

library_g4.o:	library.C nbody_particle.h
	$(CXX) -DHARP3 $(CXXFLAGS) -o library_g4.o -c library.C

nbody.o : nbody.C  nbody_particle.h design.c
	$(CXX)  -c    $(CXXFLAGS)  nbody.C 
nbody : nbody.C  nbody_particle.h $(LIBOBJS) design.c
	$(CXX) -DEVOLVE    $(CXXFLAGS) -o nbody nbody.C $(LIBOBJS) $(IONEMO) $(PGLIB) -lm

# nbody_g5 : nbody.C  nbody_particle.h $(LIBOBJSG5)
#	$(CXX)  -p -DEVOLVE -DGRAPE5 $(CXXFLAGS) $(G5INCLUDE) -o nbody_g5 nbody.C\
#	$(LIBOBJSG5) $(IONEMO) $(PGLIB) $(G5LIB) -lm

nbody_g5 : nbody.C  nbody_particle.h $(LIBOBJSG5)
	$(CXX) -DEVOLVE -DGRAPE5 $(CXXFLAGS) $(G5INCLUDE) -o nbody_g5 nbody.C\
	$(LIBOBJSG5) $(IONEMO) $(PGLIB) $(G5LIB) -lm

nbody_g4 : nbody.C  nbody_particle.h $(LIBOBJSG4)
	$(CXX) -DEVOLVE $(CXXFLAGS) -o nbody_g4 nbody.C $(LIBOBJSG4) $(IONEMO) $(PGLIB) $(G4LIB) -lm

nbody_h3 : nbody.C  nbody_particle.h $(LIBOBJSG4)
	$(CXX) -DEVOLVE    $(CXXFLAGS) -o nbody_h3 nbody.C $(LIBOBJSG4) $(IONEMO) $(H3LIB) -lm

BHtree.o : BHtree.C  nbody_particle.h BHtree.h 
	$(CXX)    $(CXXFLAGS) -c BHtree.C 
BHtree_g4.o : BHtree.C  nbody_particle.h BHtree.h 
	$(CXX) -DHARP3   $(CXXFLAGS) -c BHtree.C
	mv BHtree.o  BHtree_g4.o 
BHtree_g5.o : BHtree.C  nbody_particle.h BHtree.h 
	$(CXX) -DGRAPE5 $(G5INCLUDE) $(CXXFLAGS) -c BHtree.C
	mv BHtree.o  BHtree_g5.o 
pgetopt.o : pgetopt.C  
	$(CXX)  $(CXXFLAGS)  -c pgetopt.C
samplelog:  nbody
	nbody -i samplein > samplelog
samplelog_g4:  nbody_g4
	nbody_g4 -i samplein > samplelog_g4
samplelog_g5:  nbody_g5
	nbody_g5 -i samplein > samplelog_g5

clean:
	rm *.o

bak:	export.tar.gz
	mv export.tar.gz nbody.tar.gz

export: export.tar.gz
	ls -al export.tar.gz

export.tar.gz: $(EXPORT_FILES)
	tar cvf export.tar  $(EXPORT_FILES)
	gzip export.tar 
