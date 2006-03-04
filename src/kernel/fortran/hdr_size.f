      PROGRAM hdr
c
c     test for configure: this program writes a single integer
c     a complementing C programs reads the size of this file
c     from which the fortran header size (UNFIO_HDR_SIZE) is
c     computed.
c
      INTEGER one
      DATA one/1/

      OPEN(9,FILE='one.dat',FORM='UNFORMATTED')
      WRITE(9) one
      CLOSE(9)

      END
