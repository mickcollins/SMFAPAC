      subroutine readmoments(natomf,m,ch,dip,coord,natom1)
c     subroutine readmoments(natomf,m,lab,ch,dip,coord,natom1)
      implicit double precision(a-h,o-z)

      dimension ch(natomf),dip(natomf,3),coord(natomf,3)
      character*2 lab(natomf)
      character*20 ca, ca1

c open and read the nb.*.0.cart file

       call filelabel(m,ca1)
       n1=index(ca1,' ')-1
       ca='nb.'//ca1(1:n1)//'.0'//'.cart'

       open(unit=1,file=ca,status='old')

c read the number of atoms in this frag
       read(1,*)
       read(1,*)natom1
c start reading the data
      do n=1,natom1
       read(1,2)lab(n)
2      format(a2)
       do k=1,3
        read(1,*)coord(n,k)
       enddo
       read(1,*)
       read(1,*)ch(n)
       read(1,*)
       do k=1,3
        read(1,*)dip(n,k)
       enddo

c skip over the rest of the data for this atom
       do k=1,34
        read(1,*)
       enddo

c convert to Bohr
      Bohr=1.88972598860
      do k=1,3
       coord(n,k)=coord(n,k)*Bohr
      enddo


c end loop over atoms
      enddo

      close(unit=1)
      return
      end
