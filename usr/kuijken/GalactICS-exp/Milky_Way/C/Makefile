#Build galaxy models
BIN=../../bin

galaxy: disk bulge halo
	$(BIN)/mergerv disk bulge halo > galaxy

potential: dbh.dat

disk: freqdbh.dat cordbh.dat dbh.dat in.disk
	$(BIN)/gendisk < in.disk > disk

bulge: dbh.dat in.bulge
	$(BIN)/genbulge < in.bulge > bulge 

halo: dbh.dat in.halo
	$(BIN)/genhalo < in.halo > halo

dbh.dat: in.dbh
	$(BIN)/dbh < in.dbh

freqdbh.dat: dbh.dat h.dat b.dat
	$(BIN)/getfreqs

cordbh.dat: dbh.dat freqdbh.dat in.diskdf
	$(BIN)/diskdf < in.diskdf

clean:
	rm -f disk bulge halo galaxy

veryclean:
	rm -f *.ps *.dat disk bulge halo galaxy
