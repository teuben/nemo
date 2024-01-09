#
#   JSPAM : Interacting galaxies modeller
#   Code site:       https://github.com/jfwallin/JSPAM
#   Described in:    http://adsabs.harvard.edu/abs/2016A%26C....16...26W
#   ASCL:            https://ascl.net/1511.002


URL1 = https://github.com/jfwallin/JSPAM
URL2 = https://ascl.net/1511.002

JAR = JSPAM/archive/build/jar/jspamarchive.jar
OUT = JSPAM/archive/output

help:
	@echo "This covers:"
	@echo "JSPAM:  $(URL2)"
	@echo "   bench1      HST_NGC6786"
	@echo "   bench2      all 62 file data.txt"

JSPAM:
	git clone $(URL1)

clean:
	rm -rf $(JAR) $(OUT)


$(JAR):
	(cd JSPAM/archive; make)

# bench1 files in:  JSPAM/archive/output/HST_NGC6786
bench1: $(JAR)
	(cd JSPAM/archive ; ./runone.sh `tail -1 data.txt`)
	./spam2snap.sh JSPAM/archive/output/HST_NGC6786
	@echo "Benchmark takes ~ 1.7 sec on a nemobench5=1000 CPU"

bench2: $(JAR)
	(cd JSPAM/archive ;./runall.sh)
	@echo "Results in $(OUT)"
	@echo "Benchmark takes ~ 76 sec on a nemobench5=1000 CPU"

out:
	ls $(OUT)

JSPAM/fortran/basic_run:
	(cd JSPAM/fortran/; make)

bench3: JSPAM/fortran/basic_run
	(cd JSPAM/fortran/; ./basic_run)
	tabtos JSPAM/fortran/a.101 bench3.dat block1=pos,vel
