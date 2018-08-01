      subroutine readsigns
      use momentheader

c     open(unit=4,file='signs.out_Lev1',status='old')
c     allocate(isign(nfrag))
c     do m=1,nfrag
c     read(4,*)isign(m)
c     enddo

c     close(unit=4)

      open(unit=1,file='OUT_L1_FINALSIGNS',status='old')
      read(1,*)
      read(1,*)nL1
      allocate(isign(nL1))
      read(1,*)
      do n=1,nL1
       read(1,*)isign(n)
      enddo
      close(unit=1)


      return
      end
