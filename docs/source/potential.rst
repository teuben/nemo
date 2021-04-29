Potential(tmp)
--------------

Here we list some of the standard potentials available in NEMO,
in a variety of units, so not always :math:`G=1`!

Recall that most NEMO program use the keywords **potname=**
for the identifying name, **potpars=** for an optional
list of parameters 
and **potfile=** for an optional text string,for example
for potentials that need some kind of text file.
The parameters listed in **potpars=** will always
have as first parameter the pattern speed in cases
where rotating potentials are used. 
A Plummer potential with mass 10 and
core radius 5 would be hence be supplied 
as: ``potname=plummer potpars=0,10,5``.  The plummer potential
ignored the **potfile=** keyword.

**plummer** Plummer potential (BT, pp.42, eq. 2.47)

.. math::

    \Phi = -  {  M  \over  {   {(r_c^2 + r^2)}^{1/2} }  }



- :math:`\Omega_p`    : Pattern Speed (always the first parameter in **potpars=**)

- :math:`M`           : Total mass (clearly G=1 here)

- :math:`r_c`         : Core radius


**potname=plummer potpars=** :math:`\Omega_p`, :math:`M`, :math:`r_c`
