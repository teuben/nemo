# Branches in NEMO

Here we record some (historic) information on branches. There are a
number of branches in NEMO, some of which can be deleted, some of
which have been abandoned. Nothing was ever documented.

The first thing each branch should do is document it's use, also
provided another easy to see point where branching occured. When the
branch is merged back, it should become part of the master.

## MATDEF3

In current NEMO large 3D cubes are not well processed via FITS, since
FITS is column major (like fortran) but C (nemo) is row major.

Experimenting with different definitions of indexing multi-dimensional
arrays:

     image.h      CubeValue(i,x,y,z)    - can be different depending on MATDEF
                  jiffy notation
                  ndarray notation      - could be column or row major

