#
#



run1:
	runlpvi run1 posvel=0.2,0.34,0.5,0.1235 prefix=log dt=0.05 tstop=25000 orbit=1 ci=1::7 nstep=20
	tabplot run1/log.orb 2 3 -0.5 0.5 -0.5 0.5 yapp=1/xs
	mkorbit run1/orb.in 0.2 0.34 0 0.5 0.1235 0 potname=log potpars=0,2,0,0.7
	orbint run1/orb.in run1/orb.out  50000 0.05 ndiag=10000
	orbplot run1/orb.out xrange=-0.5:0.5 yrange=-0.5:0.5 yapp=2/xs

orb:
	rm run1/orb.*
	mkorbit run1/orb.in 0.2 0.34 0 0.5 0.1235 0 potname=log potpars=0,2,0,0.7
	orbint run1/orb.in run1/orb.out  500 0.05 ndiag=100
	orbplot run1/orb.out xrange=-0.5:0.5 yrange=-0.5:0.5 yapp=2/xs
