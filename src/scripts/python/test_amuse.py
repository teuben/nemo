#! /usr/bin/env python
#
import os, sys


try:
    import amuse
    #from amuse.community.bhtree.interface import BHTree
    from amuse.datamodel import Particles
    from amuse import datamodel
    from amuse.units import nbody_system
    from amuse.units import units
    from amuse.ic.plummer import new_plummer_sphere

    convert_nbody = nbody_system.nbody_to_si(100.0 | units.MSun, 1 | units.parsec)
    n=1000
    plummer = new_plummer_sphere(n, convert_nbody)
    stars = plummer.copy()
    print('Plummer n=',n)
    print('KE=',plummer.kinetic_energy())
    print('PE=',plummer.potential_energy())
    
except:
    print("Failing to load amuse modules")
    sys.exit(1)
