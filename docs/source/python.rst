Python
------

Python support in NEMO itself is still rudimentary, though the common

.. code-block::

     cd $NEMO
     pip3 install -e .

will add the ``nemopy`` module to your python environment. Currently only
a simple getparam interface is available.


There are good import and export routines to other N-body formats, through
which other toolkits can provide interesting ways to analyze NEMO
data. We list a few

* yt toolkit: https://yt-project.org/

* nobodykit: https://nbodykit.readthedocs.io/

* clustep, galstep: https://github.com/elvismello   create initial conditions in gadget2 format. ($NEMO/usr/mello)

* galanyl:  https://hg.sr.ht/~ngoldbaum/galaxy_analysis  pip install galanyl

* galpy:  https://github.com/jobovy/galpy

* unsio: pip install python-unsio python-unsiotools ($NEMO/usr/jcl)

* gala:  https://github.com/adrn/gala

* agama: https://github.com/GalacticDynamics-Oxford/Agama - pip install agama ($NEMO/usr/agama)

* pynbody:  https://github.com/pynbody/pynbody

* martini:  https://github.com/kyleaoman/martini ($NEMO/usr/martini)

* galpak: http://galpak3d.univ-lyon1.fr/index.html - pip install galpak

* amuse: https://github.com/amusecode   - pip install amuse-framework  

* pygad:    https://bitbucket.org/broett/pygad

* SnapGadget: https://github.com/regmachado/SnapGadget

