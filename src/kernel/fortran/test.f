C	
C	Test program for NEMO's footran interface
C		25-jun-91  1.0
C		24-may-92  1.1
C		23-apr-93  ?.?
C
C:  	Test program for NEMO's footran interface
C+
C   in=???\n		Required (dummy) filename
C   n=1000\n            Test integer value
C   pi=3.1415\n         Test real value
C   e=2.3\n             Another test value
C   text=hello world\n  Test string
C   VERSION=1.1\n       24-may-92 PJT
C-
C
      SUBROUTINE nemomain
C
C   
      include 'getparam.inc'

      INTEGER n
      DOUBLE PRECISION pi,e
      CHARACTER text*40, file*80

      file = getparam('in')
      n = getiparam('n')
      pi = getdparam('pi')
      e = getdparam('e')
      text = getparam('text')

      WRITE (*,*) 'TEST: n=',n,' pi=',pi,' e=',e,' text='//text

      RETURN
      END
