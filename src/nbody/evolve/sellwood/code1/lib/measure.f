      subroutine measure( done )
c written by Jerry Sellwood for the galaxy simulation code
c
c Rescales the global integrals and reports them out
c
c The energy, momentum, etc. (in internal program units) can be summed over all
c   particles when requested as they are stepped forward in subroutine step.
c   This routine rescales those results to external units and prints them out.
c
c calling argument
      logical done
c
c common blocks
c
      include 'admin.h'
c
c local variables
      integer i
      real pmfac, t
c
      if( done )then
        if( ncheck .eq. 0 )call
     +          crash( 'MEASURE', 'No particles in integrals check' )
c scale variables to external units - pmfac is the mass of 1 particle
        pmfac = pmass / ( lscale**3 * ts**2 )
        do i = 1, 3
          amom( i ) = amom( i ) * pmfac / ( lscale**2 * ts )
          com( i ) = com( i ) / ( lscale * real( ncheck ) )
          lmom( i ) = lmom( i ) * pmfac / ( lscale * ts )
        end do
        claus = claus * pmfac / ( lscale * ts )**2
        ke = .5 * ke * pmfac / ( lscale * ts )**2
        pe = pe * pmfac / ( lscale * ts )**2
c report
        t = ts * real( istep )
        write( no, * )
        write( no, '( '' Global integrals at time'', f12.2 )' )t
        write( no, '( 5x, ''Kinetic, potential and total E'', 3f10.4 )'
     +              )ke, pe, pe + ke
        write( no, '( 5x, ''Virial of Clausius'', f10.4 )' )claus
        write( no, '( 7x, ''Centre of mass'', 3f13.7 )' )com
        write( no, '( 6x, ''Linear momentum'', 3f13.7 )' )lmom
        write( no, '( 5x, ''Angular momentum'', 3f13.7 )' )amom
        write( no, '( 5x, ''based on'' i12, '' particles'' )' )ncheck
      else
c initialize
        do i = 1, 3
          amom( i ) = 0
          com( i ) = 0
          lmom( i ) = 0
        end do
        claus = 0
        ke = 0
        pe = 0
        ncheck = 0
      end if
      return
      end
