      subroutine readsigns
      use derivheader

      open(unit=4,file='signs.out',status='old')
      allocate(isign(nfrag))
      do m=1,nfrag
      read(4,*)isign(m)
      enddo

      close(unit=4)

      return
      end
