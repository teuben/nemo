Bonsai
------

Bonsai needs a GPU, see Makefile for example how to compile.

Units are virial units, so a NEMO plummer sphere can be integrated
without scaling the model, but it will need a conversion to/from tipsy format. Here is an example.

      mkplummer - 10000 | snaptipsy - p10k.in
      /usr/bin/time bonsai2 --infile p10k.in --snapname p10k --dev 0 --snapiter 1 -T 100 -t 0.015625 
	  # 16.97user 4.90system 0:21.78elapsed 100%CPU 
      # convert back to snapshot
      rm -f p10k.out
      for f in p10k_*.*; do
         echo $f
         tipsysnap $f - | csf - - item=SnapShot >> p10k.out
      done
      # check lagrangian radii (noting column 1 has a bug, so use the row number of the table)
      snapmradii p10k.out | tabplot - 0 2:10 line=1,1
	  
compare this with the same run using gyrfalcON:

      mkplummer p10k.snap 10000
      /usr/bin/time gyrfalcON p10k.snap p10k.gf.out tstop=100 kmax=6 step=1 theta=0.75 eps=0.05 
	  # 94.51user 0.15system 1:34.90elapsed 99%CPU
	  snapmradii p10k.gf.out | tabplot - 0 2:10 line=1,1

The real power comes from the scaling.   gyrfalcON scales O(N), so one would expect on this
computer the simulation of a 100k Plummer sphere to take 940sec.

bonsai2 took 22sec for 10k, and a 100k simulation took 33sec!

Examples taken from a 2022 desktop: Xeon(R) CPU E5-2687W 0 @ 3.10GHz;  Nvidia RTX 3070 w/ 8GB memory


