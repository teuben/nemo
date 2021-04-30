#

nbody=10000000

echo mkplummer . $nbody > run1
cat run1 run1 > run2
cat run2 run2 > run4
cat run4 run4 > run8
cat run8 run8 > run16
cat run16 run16 > run32
cat run32 run32 > run64
cat run64 run64 > run128
cat run128 run128 > run256

nemobench mkplummer . $nbody 
nemobench parallel --jobs 1 < run1
nemobench parallel --jobs 2 < run2
nemobench parallel --jobs 4 < run4
nemobench parallel --jobs 8 < run8
nemobench parallel --jobs 16 < run16
