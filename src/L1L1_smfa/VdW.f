      subroutine VdW(ndim,na,lab,radi,sym,numb)
      implicit double precision(a-h,o-z)
      dimension radi(ndim)
c     dimension numb(91)
      dimension numb(110)
c     character*2 lab(ndim),sym(91)
      character*2 lab(ndim),sym(110)

      Bohr=1.d0/1.88972598860

c get the VdW radius for each atom in Angstrom
      do n=1,na
c      do m=1,91
       do m=1,110
        if(sym(m).eq.lab(n))then
         number=m
c        radi(n)=anum2vdw(numb(m))*Bohr
         radi(n)=anum2vdw(number)*Bohr
         go to 1
        endif
       enddo
       write(*,*)' subroutine covalent failed to match'
1      continue
      enddo

      return
      end
