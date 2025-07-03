#
#
#

URL1 = https://github.com/GalacticDynamics-Oxford/Agama

agama:
	git clone $(URL1)

install: agama
	pip install -e Agama
	@echo Agama/agama.so copied to $(NEMOOBJ)/acc 

#

clean:
	rm -f sat.nemo sat_shift.nemo sat_out1.nemo sat_out2.nemo sat_out3.nemo  center.txt accel.txt

sat.nemo sat_shift.nemo:
	mkplum out=sat.nemo nbody=10000 r_s=0.1 mass=0.1 seed=123
	snapshift sat.nemo sat_shift.nemo rshift=-1,0,0 vshift=0,-0.6,0

# creates center.txt accel.txt
center.txt accel.txt:
	./mk_center.py


test1:  sat_shift.nemo
	gyrfalcON sat_shift.nemo sat_out1.nemo eps=0.02 kmax=5 step=0.25 tstop=50 accname=agama accfile=pot1.ini
	snaptrim sat_out1.nemo - times=50 | bsf - "-0.115841 1.03215 -33.2085 50 70001"



test2: center.txt sat.nemo
	gyrfalcON sat.nemo sat_out2.nemo eps=0.02 kmax=5 step=0.25 tstop=50 accname=agama accfile=pot2.ini
	snaptrim sat_out2.nemo - times=50 | bsf - "-4.50823 10.815 -36.6663 50 70001"


#   convert back to host centered coordinate system
test2b:
	snapprint sat_out2.nemo x,y,z header=time |\
	awk '{print $$2-cos(0.6*$$1),$$3-sin(0.6*$$1),$$4;}' |\
	tabtos - sat_out3.nemo nbody=10000 header=time block1=x,y,z


# see also more elaborate code in $AGAMA_DIR/py/example_nbody_simulation.py
