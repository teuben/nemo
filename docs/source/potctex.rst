Example Potentials
------------------



In the latex manual this chapter is derived from the code. We thus need a new ``crst`` script
that produces this.   Needs to handled embedded math, such as :math:`\frac{ \sum_{t=0}^{N}f(t,k) }{N}`


bar83
~~~~~
% File: bar83.c

potname=bar83 potpars={\it $\Omega,f_m,f_x,{c\over a}$}} 

Barred potential as described by Teuben and Sanders (1983), 
see also {\bf teusan83}. 

Note: the potential only valid in the z=0 plane! 


bulge1
~~~~~~
% File: bulge1.c

{\bf potname=bulge1 
potpars={\it $\Omega,M,R,c/a$}} 

homogeneous oblate bulge with mass $M$, radius $R$, and axis ratio $c/a$ 

ccd
~~~
% File: ccd.c

{\bf potname=ccd 
potpars={\it $\Omega,Iscale,Xcen,Ycen,Dx,Dy$} 
potfile={\it image(5NEMO)}} 

This potential is defined using a simple cartesian grid on which the potential 
values are stored. Using bilinear interpolation the values and derivatives are 
computed at any point inside the grid. Outside the grid (as defined by the 
WCS in the header) the potential is not defined and assumed 0. 
The lower left pixel of an image in NEMO is defined as (0,0), with WCS values 
Xmin,Ymin derived from the header. If the (Xcen,Ycen) parameters are used, 
these are the 0-based pixel coordinates of the center pixel. If (Dx,Dy) are used, 
these are the pixel separations. 
To aid astronomical images where $Dx < 0$, these are interpreted as positive. 
Also note that potentials are generally negative, so it is not uncommon to need 
$Iscale = -1$. Programs such as {\it potccd} can create such a {\bf ccd} grid 
potential from a regular potential. 

Note: Since these forces are defined only in the Z=0 plane, the Z-forces are always 
returned as 0. 

cp80
~~~~
% File: cp80.c

{\bf potname=cp80 
potpars={\it $\Omega,\epsilon$}} 

Contopoulos \& Papayannopoulos (1980, A\&A, 92,33) 
used this potential 
in the study of orbits in barred galaxies. Note that their 
*bar* is oriented along the Y-axis, an axis ratio is not 
well defined, and for larger values of $\epsilon$ the density 
can be negative. The potential used is given by adding an 
axisymmetric component to a m=2 fourier component: 

.. math::

   \Phi = \Phi_1 + \Phi_2 

where :math:`\Phi_1` is the Isochrone potential with unit scalelength and 
mass, and :math:`\Phi_2` the Barbanis & Woltjer (1965) potential:

.. math::

   \Phi_1 = - { 1 \over { (1 + \sqrt{1+r^2})}} 

and

.. math::

   \Phi_2 = \epsilon r (16-r) cos(2\phi) 


A value of :math:`\epsilon=0.00001` is the default for a moderate bar, 
whereas 0.001 is a strong bar! 

dehnen
~~~~~~
% File: dehnen.c

{\bf potname=dehnen 
potpars={\it $\Omega,M,a,\gamma$}} 

Walter Dehnen (1993, MN {\bf 265}, 250-256) introduced a 
family of potential-density pairs for spherical systems: 

The potential is given by: 

.. math::

   \Phi = { G M \over a } {1\over{2-\gamma}} {\left[ 1 - {\left(r\over{r+a}\right)}^{2-\gamma}\right]} 

cumulative mass by 

.. math::

   M(r) = M { r \over {(r+a)}^{3-\gamma} } 

and density by 

.. math::

   \rho = { {(3-\gamma)M} \over {4\pi}} { a \over {r^{\gamma} (r+a)^{4-\gamma}}} 

with $0 <= \gamma < 3$. 
Special cases are the Hernquist potential ($\gamma=1$), and the 
Jaffe model ($\gamma=2$). The model with $\gamma=3/2$ seems to 
give the best comparison withe de Vaucouleurs $R^{1/4}$ law. 

See also Tremaine et al. (1994, AJ, 107, 634) in which they describe 
the same density models with $\eta=3-\gamma$ and call them 
$\eta$-models. 

dublinz
~~~~~~~
% File: dublinz.c

{\bf potname=dublinz 
potpars={\it $\Omega,r_0,r_1,v_1,dvdr,s,h$}} 

Forces defined by a double linear rotation curve defined by 
($r_1,v_1$) and a gradient $dvdr$ between $r_0$ and $r_1$. 
As in {\bf flatz} (from which this one is derived), the 
potential is quasi harmonic in $Z$ (linear forces), 
with radial scalelength $h$ and scale height $s$ 

expdisk
~~~~~~~
% File: expdisk.c

{\bf potname=expdisk 
potpars={\it $\Omega,M,a$}} 

Exponential disk (BT, pp.77)

.. math::

   \Phi = - {M \over r_d} x \left[ I_0(x)K_1(x) - I_1(x)K_0(x) \right] 


flatz
~~~~~
% File: flatz.c

potname=flatz potpars=:math:`\Omega,r_0,v_0,s,h`

forces defined by a rotation curve that is linear to 
$(r_0,v_0)$ and flat thereafter and quasi harmonic in $Z$, 
with radial scalelength $h$ and scale height $s$. 
See also {\bf dublinz} for a variation on this theme. 


halo
~~~~
% File: halo.c

{\bf potname=halo 
potpars={\it $\Omega,v_0,r_c$}} 


hh64
~~~~
% File: hh64.c

potname=hh64   potpars=:math:`\Omega,\lambda`

.. math::

       \Phi = {1 \over 2} ( x^2 + x^2 ) + \lambda ( x^2 y - {1\over 3} y^3 )



grow_plum
~~~~~~~~~
% File: grow_plum.c


grow_plum2
~~~~~~~~~~
% File: grow_plum2.c


harmonic
~~~~~~~~
% File: harmonic.c

{\bf potname=harmonic 
potpars={\it $\Omega,\omega_x^2,\omega_z^2,\omega_z^2$}} 


Harmonic potential 

.. math::

    \Phi = {1 \over 2} \omega_x^2 x^2 + {1 \over 2} \omega_y^2 y^2 + {1 \over 2} \omega_z^2 z^2 


hernquist
~~~~~~~~~
% File: hernquist.c

{\bf potname=hernquist 
potpars={\it $\Omega,M,r_c$}} 

The Hernquist potential (ApJ, 356, pp.359, 1990) is a special $\gamma=1$ case 
of the Dehnen potential. The potential is given by:

.. math::

   \Phi = - { M \over {(r_c+r)}} 

and mass 

.. math::

   M(r) = M { r^2 \over {(r+r_c)}^2 } 

and density 

.. math::

   \rho = { M \over {2\pi}} {r_c \over r} { 1 \over {(r+r_c)}^3} 


hom
~~~
% File: hom.c

{\bf potname=hom  potpars={\it $\Omega,M,R,\tau$}} 

hubble
~~~~~~
% File: hubble.c

{\bf potname=hubble 
potpars={\it $\Omega,M,R,b,c$}} 
where $M$ and $R$ are the core mass and radius. $b$ and $c$ are, if 
given, the intermediate and short axes can be different from the 
core radius. 

The Hubble profile (BT, pp 39, req. 2-37 and 2-41) has a density 
law: 

.. math::

   \rho = \rho_h ( 1 + (r/r_h)^2 )^{-3/2} 

and an equally simple expression for the projected surface brightness: 

.. math::

   \Sigma = 2 \rho_h r_h ( 1 + (r/r_h)^2)^{-1} 

The derivation of the potential is a bit more involved, since there 
is no direct inversion, and integration in parts is needed. The 
cumulative mass is given by:

.. math::

   M_h(r) = 4\pi r_h^3 \rho_h \{ \ln[(r/r_h) + \sqrt{1+(r/r_h)^2}] - { {r/a} \over { \sqrt{1+(r/r_h)^2} } } \} 

and the potential

.. math::

   \Phi(r) = - { {GM_h(r)}\over {r} } - { {4\pi G \rho_h r_h^2} \over {\sqrt{1+r}} } 


kuzmindisk
~~~~~~~~~~
% File: kuzmindisk.c

{\bf potname=kuzmin 
potpars={\it $\Omega,M,a$}} 

Kuzmin (1956) found a closed expression for the potential of 
an infinitesimally thin disk with a Plummer potential in the 
plane of the disk (see also BT pp43, eq. 2-49a and 2-49b): 

.. math::

   \Phi = - { G M \over {\sqrt{r^2 + (a+{|z|})^2}}} 

and corresponding surface brightness ({\it check units}) 


.. math::

   \Sigma = { {a M} \over {2 \pi {(a^2 + r^2)}^{-3/2}}} 

With $GMa^2 = V_0^2$. 
This potential is also known as a Toomre n=1 disk, since it 
was re-derived by Toomre (1963) as part of a series of disks 
with index $n$, where this disk has $n=1$. 

isochrone
~~~~~~~~~
% File: isochrone.c

{\bf potname=isochrone 
potpars={\it $\Omega,M,R$}} 

jaffe
~~~~~
% File: jaffe.c

{\bf potname=jaffe 
potpars={\it $\Omega,M,r_c$}} 

The Jaffe potential (BT, pp.237, see also MNRAS 202, 995 (1983))), 
is another special $\gamma=2$ case of the Dehnen potential. 

.. math::

   \Phi = - { M \over r_c} \ln{ \left( { r \over {r_c + r} } \right) } 


log
~~~
% File: log.c

% CTEX Line: 8
{\bf potname=log 
potpars={\it $\Omega,M_c,r_c,q$}} 

The Logarithmic Potential (BT, pp.45, eq. 2.54 and eq. 3.77) has 
been often used in orbit calculations because of its flat rotation 
curve. The potential is given by 

.. math::

   \Phi = {1\over 2} v_0^2  \ln{ \left( r_c^2 + r^2 \right) } 


with $ M_c \equiv {1\over 2} r_c v_0^2 $ defined as the *core mass*. 

mestel
~~~~~~

% File: mestel.c

% CTEX Line: 10
{\bf potname=mestel 
potpars={\it $\Omega,M,R$}} 

miyamoto
~~~~~~~~
% File: miyamoto.c

% CTEX Line: 20
{\bf potname=miyamoto 
potpars={\it $\Omega,a,b,M$}} 

.. math::


   \Phi = - { M \over { .... } }



nfw
~~~
% File: nfw.c
% CTEX Line: 29

The NFW (Navarro,Frank \& White) density is given by 

.. math::

   \rho = { M_0 \over { r (r+a)^2}} 


and the potential by

.. math::

   \Phi = -4 \pi M_0 { \ln{(1+r/a)} \over r } 


null
~~~~
% File: null.c

% CTEX Line: 5

This potential has no other meaning other than to fool the compiler. 
It has no associates potential, thus the usual potname, potpars,potfile 
will have no meaning. Use {\bf potname=zero} if you want a real potential 
with zero values. 

op73
~~~~
% File: op73.c

% CTEX Line: 14
{\bf potname=op73 
potpars={\it $\Omega,M_H,r_c,r_h$}} 

Ostriker-Peebles 1973 potential 
(1973, ApJ {\bf 186}, 467). 
Their potential is given in the form of the radial force law in the disk 
plane: 

.. math::

   F = { M \over R_h^2 } { {(R_h+R_c)}^2 \over {(r+R_c)}^2 } { r \over R_h } 


plummer
~~~~~~~
% File: plummer.c

% CTEX Line: 8
{\bf potname=plummer 
potpars={\it $\Omega,M,R$}} 

Plummer potential (BT, pp.42, eq. 2.47, see also MNRAS 71, 460 (1911)) 

.. math::

   \Phi = - { M \over { {(r_c^2 + r^2)}^{1/2} } } 


plummer2
~~~~~~~~
% File: plummer2.c

rh84
~~~~
% File: rh84.c

% CTEX Line: 20
{\bf potname=rh84 
potpars={\it $\Omega,B,a,A,r_0,i_0,j$}} 


This 2D spiral and bar potential was used by Robert and collaborators 
in the 70s and 80s. 
For counterclockwise streaming, this spiral is a trailing 
spiral when the pitch angle ($i_0$) is positive. 
Within a radius $r_0$ the potential becomes barlike, with 
the bar along the X axis. 
At large radii the spiral is logarithmic. 
References: 

Roberts \& Haussman (1984: ApJ 277, 744) 

Roberts, Huntley \& v.Albada (1979: ApJ 233, 67) 

rotcur0
~~~~~~~
% File: rotcur0.c

% CTEX Line: 9
{\bf potname=rotcur0 
potpars={\it $\Omega,r_0,v_0$}} 

The forces returned are the axisymmetric forces as defined by 
a linear-flat rotation curve as defined by the turnover point $r_0,v_0$. 
The potential is not computed, instead the interpolated rotation 
curve is returned in as the potential value. 

rotcur
~~~~~~
% File: rotcur.c

% CTEX Line: 14
{\bf potname=rotcur 
potpars={\it $\Omega$} 
potfile={\it table(5NEMO)}} 

The forces returned are the axisymmetric forces as defined by 
a rotation curve as defined by a table given from an ascii table. 
The potential is not computed, instead the interpolated rotation 
curve is returned in as the potential value. 

This version can only compute one version; i.e. 
on re-entry of inipotential(), old versions are lost. 

sh76
~~~~
% File: sh76.c

{\bf potname=sh76
       potpars={\it $\Omega,A,\alpha,\epsilon$}}

This bar potential was used by Sanders and Huntley (1976) and
also used in Sanders (2019).   The density perturbation is given
by

.. math::

   \sigma(r,\theta) = A r^{-\alpha} (1+\epsilon\cos{2\theta})

and the potential

.. math::

    \Phi(r,\theta) = -2\pi G c_1 A r^{-\alpha+1} {1 \over {1-\alpha}} ( 1 + \beta (\alpha-1) \cos{2\theta})

where

.. math::

   \beta =  { {(2-\alpha)} \over { \alpha(3-\alpha)} }  \epsilon

.. math::

   c_1 = { { \Gamma{[{1\over 2}(2-\alpha)]}  \Gamma{[{1\over 2}(\alpha+1)]} }   \over
          { \Gamma{[{1\over 2}\alpha]}  \Gamma{[{1\over 2}(3-\alpha)]} } }


teusan85
~~~~~~~~
% File: teusan85.c

% CTEX Line: 25
{\bf potname=teusan85} 

This potential is that of a barred galaxy model as 
described in Teuben \& Sanders (1985) 
This bar is oriented along the X axis. 
This is the 2D version for forces. This version should give (near) 
identical results to {\bf bar83} and very simlar to {\bf athan92}. 





triax
~~~~~
% File: triax.c

% CTEX Line: 11
{\bf potname=triax} 

A growing bi/triaxial potential 



twofixed
~~~~~~~~
% File: twofixed.c

% CTEX Line: 16
{\bf potname=twofixed 
potpars={\it $\Omega,M_1,x_1,y_1,z_1,M_2,x_2,y_2,z_2$}} 


This potential is defined by two fixed points, with different masses 
and positions. Orbits in this potential exhibit a number of interesting 
properties. One well known limit is the {\tt stark problem}, where one 
of the two bodies is far from the other and near-circular orbits near 
the central particles are studied. Another is the limit or two particles 
near to other and orbits that circumscribe both particles. 



plummer4
~~~~~~~~
% File: plummer4.c

% CTEX Line: 10
potname=plummer potpars=:math:`\Omega,M,R`

Plummer potential (BT, pp.42, eq. 2.47, see also MNRAS 71, 460 (1911)) 

.. math::

   \Phi = - { M \over { {(r_c^2 + r^2)}^{1/2} } } 


vertdisk
~~~~~~~~
% File: vertdisk.c


tidaldisk
~~~~~~~~~
% File: tidaldisk.c
% CTEX Line: 8

Tidal field exerted by a (plane-parallel) stellar disk on a cluster passing 
through with constant vertical velocity. 
Useful for simulations of disk-shocking of, say, globular clusters 

The following three density models are available 

1. thin disk: 

.. math::

   \rho(z) = \Sigma \delta(z) 

2. exponential disk: 

.. math::

   \rho(z) = {\Sigma \over {2h}} \exp{ { -|z|} \over h} 


3. sech$^2$ disk: 

.. math::

   \rho(z) = {\Sigma \over {4h}} sech^2{ { z \over {2h}}} 


Parameters (to be given by potpars=...) are: 

.. sourcecode::

    par[0] = not used (reserved for pattern speed in NEMO) 
    par[1] = h scale-height par[1] = 0 -> thin disk 
    par[1] > 0 -> vertically exponential disk 
    par[1] < 0 -> sech^2 disk with h=|par[1]| 
    par[2] = Sig disk surface density 
    par[3] = Vz constant vertical velocity of cluster center 
    par[4] = Z0 cluster center z-position at t=0 
    par[5] = add boolean: add tidal potential or not? 


We always assume G=1. 

If you want to include the acceleration of the disk on the cluster as a 
whole, rather than assume a constant velocity, use vertdisk.c 

Some words on the mechanics 

Assume that the plane-parallel disk potential and force are given by 

.. math::

   \Phi(Z) , F(Z) = -\Phi'(Z). 

Then, the tidal force exerted on a star at position z w.r.t. to cluster 
center, which in turn is at absolute height Zc = Z0 + t Vz, is simply 

.. math::

   F_t(z) = F(Zc+z) - F(Zc). 

Integrating this from z=0 to z gives the associated tidal potential as 

.. math::

   \Phi_t(z) = \Phi(Zc+z) - \Phi(Zc) + z  F(Zc). 

Whenever the tidal force \& potential are desired at a new time t, we 
pre-compute $Zc$ and the plane-parallel potential and force at $Z=Zc$. 
Note that when both $Zc$ and $Zc+z$ are outside of the mass of the disk (and 
$Z=0$ 
is not between them), both tidal force and potential vanish identically. 
The standard treatment of tidal forces corresponds to approximating (2) by 
$F(Zc) + z * F'(Zc)$. This method, however, breaks down for disks that are 
thin compared to the cluster, while our method is always valid, even for a 
razor thin disk. 

polynomial
~~~~~~~~~~
% File: polynomial.c

% CTEX Line: 9
{\bf potname=polynomial 
potpars={\it $\Omega,a0,a1,a2,a3,....$}} 

Polynomial potential 

.. math::

   \Phi = a_0 + a_1 r + a_2 r^2 + .... a_N r^N 


where any unused coefficients will be set to 0. Up to 16 (defined 
as MAXPOW) can be used. 

wada94
~~~~~~
% File: wada94.c

% CTEX Line: 11
{\bf potname=wada94 
potpars={\it $\Omega,c,a,\epsilon$}} 

Wada (1994, PASJ 46, 165) and also 
Wada \& Have (1992, MN 258, 82) 
used this potential in the study of gaseous orbits in barred galaxies. 

.. math::

   \Phi = \Phi_0 + \Phi_b 

where $\Phi_1$ is the Toomre potential with scalelength $a$

.. math::

   \Phi_0 = - { 1 \over \sqrt{R^2 + a^2}} 

and

.. math::

   \Phi_b = -\epsilon { {a R^2} \over { {(R^2 + a^2)}^2 } } 

A relationship for the axisymmetric component is

.. math::
   
-\sqrt(27/4) 


zero
~~~~
% File: zero.c

% CTEX Line: 6
{\bf potname=zero} 

Zero potential 

.. math::

   \Phi = 0 


