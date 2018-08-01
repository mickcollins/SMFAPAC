      program distinctlabels
      implicit double precision(a-h,o-z)

      character*2, allocatable :: lab(:)

      character*2, atom

c reads in a molecule coordinate file and works out
c the unique element labels

c the coord file is read from standard input

      read(5,*)natom

      allocate(lab(natom))

      read(5,*)

      ic=0
      do n=1,natom
       read(5,*)atom
       if(ic.gt.0)then
        do m=1,ic
         if(atom.eq.lab(m))go to 1
        enddo
c new type
        ic=ic+1
        lab(ic)=atom
       endif
1       continue
       if(ic.eq.0)then
        ic=ic+1
        lab(ic)=atom
       endif
      enddo

      open(unit=1,file='labels',status='unknown')
      do n=1,ic
       write(1,100)lab(n)
      enddo
100   format(a2)
      close(unit=1)
      end
