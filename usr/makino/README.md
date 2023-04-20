

FDPS -    framework for Developing Particle Simulator

See also https://github.com/FDPS/FDPS/

# command line flags for nbody.out:

     o: dir name of output (default: ./result)
     t: theta (default: 0.5)
     T: time_end (default: 10.0)
     s: time_step (default: 1.0 / 128.0)
     d: dt_diag (default: 1.0 / 8.0)
     D: dt_snap (default: 1.0)
     l: n_leaf_limit (default: 8)
     n: n_group_limit (default: 64)
     N: n_tot (default: 1024)
     h: help


# bench:  nbody.out -N 8192 -n 256

test1:  58.12user 0.07system 0:58.24elapsed  99%CPU 
test2: 157.69user 0.33system 0:20.64elapsed 765%CPU 
test3:  dnc
test4:  12.58user 0.12system 0:12.75elapsed 99%CPU 
test5:  dnc

(dnc = did not compile)

# comparing bench w/ gyrfalcON ?

How does gyrfalcON compare to nbody.out

    mkplummer p8k 8192
    gyrfalcON p8k p8k.out tstop=10 step=1 theta=0.5 kmax=7 eps=0.05
    19.67user 0.00system 0:19.69elapsed 99%CPU 

# input files for nbody.out ?

Alas, there is no input option for nbody.out (hence the name?), neither is there an
option for softening, but once there is, here's an example how to create that
format with NEMO

n=1000
t=0.0
nemoinp $t,$n > data.in
mkplummer - $n | snapprint - i,m,x,y,z,vx,vy,vz >> data.in
