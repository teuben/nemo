
all:	test1 test2 test3 test4 test5 test6

#  standard benchmarks
clean:
	rm -f bench1?.* 

bench1a:
	time ./plummer_stats.csh tmp=bench1a nexp=100 nbody=0 mlo=1 mhi=1 quick=1
	@echo "R_v above should show ...0.966369 0.0685018"
	@echo "tail should show ..."
	sort -n bench1a.100.tab | tail -2
	@echo Should be:
	@echo 2.51094 0.01 -0.0920454 0.247735 0.0795859 0.264282 
	@echo 3.67919 0.01 -0.18561 0.123529 -0.143081 0.222958 

bench1b:
	time ./cluster_stats.py  out='"bench1b"' nexp=100 nbody=0 mlo=1 mhi=1 quick=1
	#@echo "R_v,R_c above should show ..."
	#@echo "r_v: 0.9663689500000001 +/- 0.06850180731183303"
	#@echo "r_c: 0.3959586 +/- 0.13051354191592532"
	sort -n bench1b.100.tab | tail -2
	@echo Should be:
	@echo 2.51094 0.01 -0.0920454 0.247735 0.0795859 0.264282 
	@echo 3.67919 0.01 -0.18561 0.123529 -0.143081 0.222958

bench1c:
	@rm -f bench1c.dat
	time mkplummer bench1c.dat nbody=100 nmodel=100 mode=0 rfrac=5.0 seed=1 'massname=n(m)' masspars=p,-2.0 massrange=1,1
	time mkplummer bench1c.dat nbody=100 nmodel=100 mode=0 rfrac=5.0 seed=1

bench1d:
	mkplummer - nbody=100 nmodel=100 mode=1 rfrac=5.0 seed=1 |\
		snapsort bench1d.dat - r2 |\
		snapprint - r2,m,vx,vy,vz,v2 debug=-1 | tail -2
	@echo Should be:
	@echo 2.51094 0.01 -0.0920454 0.247735 0.0795859 0.264282 
	@echo 3.67919 0.01 -0.18561 0.123529 -0.143081 0.222958

bench1c-old:
	time ./cluster_stats.py  out='"bench1c"' nexp=100 mlo=1 mhi=1 mode=1
	time ./cluster_stats.py  out='"bench1c"' nexp=100 mlo=1 mhi=1 mode=2

bench1d-old:
	time ./cluster_stats.py  out='"bench1d"' nexp=2000 mlo=1 mhi=1 mode=1 


# below here private tests that should probably be removed

# NEMO's mkplummer
test1:
	plummer_stats.csh tmp=exp1 nexp=100 mlo=1 mhi=1
	plummer_stats.csh tmp=exp2 nexp=100 mlo=1 mhi=1 tstop=10
	plummer_stats.csh tmp=exp3 nexp=100 
	plummer_stats.csh tmp=exp4 nexp=100 tstop=10
	plummer_stats.csh tmp=exp5 nexp=100 nsmall=1000
	plummer_stats.csh tmp=exp6 nexp=100 nsmall=1000 rcut=20

# mcluster's plummer
test2:
	plummer2_stats.csh tmp=exp11 nexp=100 mlo=1 mhi=1
	plummer2_stats.csh tmp=exp12 nexp=100
	plummer2_stats.csh tmp=exp13 nexp=100 nsmall=0 msmall=232
	plummer2_stats.csh tmp=exp14 nexp=100 nsmall=0 msmall=300 rc=0.25 rbin=0:2.5:0.125
	plummer2_stats.csh tmp=exp15 nexp=100 nsmall=0 msmall=300 rc=0.25 rbin=0:2.5:0.125 shift=0.5
	plummer2_stats.csh tmp=exp16 nexp=100 nsmall=0 msmall=1000 rbin=0:4:0.25 mlo=1 mhi=1

# mcluster's homsph
test3:
	plummer3_stats.csh tmp=exp20 nexp=100 nsmall=0 msmall=1000 rc=2.5 rbin=0:2.5:0.125 mlo=1 mhi=1 
	plummer3_stats.csh tmp=exp21 nexp=100 nsmall=0 msmall=1000 rc=2.5 rbin=0:2.5:0.125 mlo=1 mhi=1 tstop=10
	plummer3_stats.csh tmp=exp22 nexp=100 nsmall=0 msmall=1000 rc=2.5 rbin=0:2.5:0.125 tstop=10
	plummer3_stats.csh tmp=exp23 nexp=100 nsmall=0 msmall=1000 rc=2.5 rbin=0:2.5:0.125 

# mcluster's king
test4:
	./plummer4_stats.csh tmp=exp30 nexp=100 nsmall=0 msmall=100 rc=0.25 rbin=0:2.5:0.125 w0=1
	./plummer4_stats.csh tmp=exp31 nexp=100 nsmall=0 msmall=100 rc=0.25 rbin=0:2.5:0.125 w0=5
	./plummer4_stats.csh tmp=exp32 nexp=100 nsmall=0 msmall=100 rc=0.25 rbin=0:2.5:0.125 w0=10 (11 ok, 12 fails)
	./plummer4_stats.csh tmp=exp33 nexp=100 nsmall=0 msmall=100 rc=0.25 rbin=0:2.5:0.125 w0=1 tstop=10
	./plummer4_stats.csh tmp=exp34 nexp=100 nsmall=0 msmall=100 rc=0.25 rbin=0:2.5:0.125 w0=10 tstop=10
	./plummer4_stats.csh tmp=exp35 nexp=100 nsmall=0 msmall=1000 rc=0.25 rbin=0:2.5:0.125 w0=1 tstop=10
	./plummer4_stats.csh tmp=exp36 nexp=100 nsmall=0 msmall=1000 rc=0.25 rbin=0:2.5:0.125 w0=10 tstop=10
	./plummer4_stats.csh tmp=exp37 nexp=100 nsmall=0 msmall=1000 rc=0.5 rbin=0:2.5:0.125 w0=1 mlo=1 mhi=1
	./plummer4_stats.csh tmp=exp38 nexp=100 nsmall=0 msmall=1000 rc=0.5 rbin=0:2.5:0.125 w0=7 mlo=1 mhi=1

# NEMO's king
test5:
	./plummer5_stats.csh tmp=exp40 nexp=100 nsmall=100 rbin=0:2.5:0.125 w0=1
	./plummer5_stats.csh tmp=exp41 nexp=100 nsmall=100 rbin=0:2.5:0.125 w0=7
	./plummer5_stats.csh tmp=exp42 nexp=100 nsmall=100 rbin=0:2.5:0.125 w0=1 tstop=10
	./plummer5_stats.csh tmp=exp43 nexp=100 nsmall=100 rbin=0:2.5:0.125 w0=7 tstop=10
	./plummer5_stats.csh tmp=exp44 nexp=100 nsmall=1000 rbin=0:2.5:0.125 w0=1 tstop=10
	./plummer5_stats.csh tmp=exp45 nexp=100 nsmall=1000 rbin=0:2.5:0.125 w0=7 tstop=10
	./plummer5_stats.csh tmp=exp46 nexp=100 nsmall=1000 rbin=0:2.5:0.125 w0=1 
	./plummer5_stats.csh tmp=exp47 nexp=100 nsmall=1000 rbin=0:2.5:0.125 w0=7 


# NEMO's mkhomsph
test6:
	./plummer6_stats.csh tmp=exp50 nexp=100 nsmall=100  rbin=0:2.5:0.125 rc=2
	./plummer6_stats.csh tmp=exp51 nexp=100 nsmall=100  rbin=0:2.5:0.125 rc=2 tstop=10
	./plummer6_stats.csh tmp=exp52 nexp=100 nsmall=1000 rbin=0:2.5:0.125 rc=2 tstop=10
	./plummer6_stats.csh tmp=exp53 nexp=100 nsmall=1000 rbin=0:2.5:0.125 rc=2
	./plummer6_stats.csh tmp=exp54 nexp=100 nsmall=1000 rbin=0:2.5:0.125 rc=2  

# NEMO's mkking (falcON)
test7:
	./plummer7_stats.csh tmp=exp60 nexp=100 nsmall=1000 rbin=0:2.5:0.125 w0=1
	./plummer7_stats.csh tmp=exp61 nexp=100 nsmall=1000 rbin=0:4.0:0.25  w0=7
	./plummer7_stats.csh tmp=exp62 nexp=100 nsmall=1000 rbin=0:20.0:0.25 w0=11
