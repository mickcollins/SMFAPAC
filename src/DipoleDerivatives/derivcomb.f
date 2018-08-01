      program combdipderiv

      use derivheader

c this program reads in the forces and (possibly) hessians
c for each of the main fragments and combines them
c into gradients and hessians for the atoms in the
c whole molecule

      real*8, allocatable  :: dvdx(:,:),fa(:,:)
      real*8 energy,diff

      mcount=0
      rcount=0.d0

      bohr=1.d0/1.8897259886d0

c get the total number of atoms in the molecule from name.xyz
      open(unit=1,file='name.xyz',status='old')
      read(1,*)natom
      close(unit=1)

c read in the atom numbering of the fragments
      call readfrag

c read in the coefficients for each fragment in the energy sum
      call readsigns

c read in the dipole derivatives for each fragment,
c and assign to the real atoms
      call readdata

c output the dipole derivative tensor
      do n=1,natom
      do k=1,3
       write(6,100)(dipdr(j,k,n),j=1,3)
      enddo
      enddo
100   format(3e14.6)

      end
