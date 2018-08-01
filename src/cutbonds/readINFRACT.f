      subroutine readINFRACT

      use fractdata
      implicit double precision(a-h,o-z)

      open(unit=1,file='IN_FRACT',status='old')
      read(1,90)icomm
      read(1,*)Level
      read(1,90)icomm
      read(1,*)dtol

c dtol is not used in the fragmentation
c but in a later program

90    format(a80)
      close(unit=1)

      return
      end
