      program combdipderiv

      use momentheader

c this program reads in the multipole moments of the atoms
c in each of the Level 1 fragments; allocates each to the
c corresponding atom in the molecule; and outputs a file
c that lists the multipole moments for each atom in the molecule
c
c These are used elsewhere to generate the ESP around the molecule.

      real*8, allocatable  :: dvdx(:,:),fa(:,:)
      real*8 energy,diff

      mcount=0
      rcount=0.d0

      bohr=1.d0/1.8897259886d0

c get the total number of atoms in the molecule from name.xyz
      open(unit=1,file='name.xyz',status='old')
      read(1,*)natom
      allocate(c(natom,3))
      allocate(lab(natom))
      allocate(ch(natom))
      allocate(d(natom,3))
      allocate(q(natom,3,3))
      allocate(o(natom,3,3,3))
      allocate(h(natom,3,3,3,3))

      do n=1,natom
       read(1,*)lab(n),(c(n,k),k=1,3)
       do k=1,3
        c(n,k)=c(n,k)/bohr
       enddo
      enddo
      close(unit=1)
c read in the atom numbering of the fragments
      call readfrag
c read in the coefficients for each fragment in the energy sum
      call readsigns

c reduce nfrag to nL1, the number of Lev1 frags after cancellation
      nfrag=nL1

c read in the dipole derivatives for each fragment,
c and assign to the real atoms
      call readdata

      open(unit=1,file='allmoments',status='unknown')
      write(1,*)'The electrostatic moments for all atoms at Level 1'
      write(1,*)natom
      do i=1,natom
       write(1,101)lab(i)
101   format(a2)
       do k=1,3
        write(1,*)c(i,k)
       enddo
       write(1,*)' Charge'
       write(1,*)ch(i)
       write(1,*)' Dipole'
       do k=1,3
        write(1,*)d(i,k)
       enddo
       write(1,*)' Quodrupole'
       do k1=1,3
        do k2=k1,3
         write(1,*)q(i,k1,k2)
         enddo
       enddo
       write(1,*)' Octapole'
        do k1=1,3
        do k2=k1,3
        do k3=k2,3
        write(1,*)o(i,k1,k2,k3)
        enddo
        enddo
        enddo
       write(1,*)' Hexadecapole'
        do k1=1,3
        do k2=k1,3
        do k3=k2,3
        do k4=k3,3
        write(1,*)h(i,k1,k2,k3,k4)
        enddo
        enddo
        enddo
        enddo
c end the i loop
      enddo
100   format(3e14.6)
      close(unit=1)

      end
