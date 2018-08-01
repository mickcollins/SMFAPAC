      subroutine readdata
      use derivheader

c the force and (possibly) hessian data for each fragment is 
c read from standard input 

      integer, allocatable   :: numberf(:)

      bohr=1.d0/1.8897259886d0

c read nflag=1 for first dervs only or nflag=2 for seconds also
      read(5,*)nflag
c     if((nflag-1).gt.1.or.(nflag-1).lt.0)then
c     write(6,*)" The flag in the input data is neither 1 or 2"
c     stop
c     endif

c get the total number of atoms in the molecule from name.xyz
      open(unit=1,file='name.xyz',status='old')
      read(1,*)natom
      close(unit=1)

c allocate some arrays
      allocate(amas(natom))
      allocate(coord(natom,3))
      allocate(lab(natom))
c amas will contain the mass of each atom
c coord will contain the Cartesian coordinates of each atom
c lab will contain the elemental symbol of each atom

c now read the data
      read(5,*)nfrag
c nfrag = the number of fragments

c allocate some arrays
      allocate(energyfrag(nfrag))
c energyfrag = the energy of each fragment
      allocate(natomfrag(nfrag))
c natomfrag = the number of atoms in a fragment (real+caps)
      allocate(nat0(nfrag))
c nat0 = the number of real atoms in the fragment
c (used in later subroutines)

      allocate(ncap(nfrag))
c ncap = the number of caps in the fragment
c (used in later subroutines)
      allocate(natcap(nfrag,maxcaps,2))
c natcap contains the (whole molecule) atom numbers of the two
c atoms whose positions define the position of the cap
      allocate(fact(nfrag,maxcaps))
c fact = the factor which allocates the fraction of a cap force
c to each of the two real atoms
      allocate(isign(nfrag))
c isign = the coefficent of the fragment in the energy sum

      allocate(numberf(nfrag))
c numberf = a temporary array for the number of atoms in a fragment

c read in the data a first time, just to get the max number
c of atoms in a fragment, so a dimension can be allocated

      do m=1,nfrag
      read(5,*)numberf(m)
c     do n=1,2*numberf(m)+5
      if(nflag.eq.0)then
       do n=1,numberf(m)+1
        read(5,*)
       enddo
      endif
      if(nflag.ge.1)then
      do n=1,2*numberf(m)+1
       read(5,*)
      enddo
      endif
      if(nflag.eq.2)then
       n3=3*numberf(m)

       do i=1,n3
       do j=1,i
        read(5,*)
       enddo
       enddo
      endif
c end the m loop
      enddo

      maxatom=0
      do m=1,nfrag
       if(numberf(m).gt.maxatom)maxatom=numberf(m)
      enddo

c now have the max number of atoms in a fragment
      m3=3*maxatom

      allocate(c(nfrag,maxatom,3))
c c = the Cartesian coords of each atom in the fragment
      allocate(force(nfrag,maxatom,3))
c force = the Cartesian forces of each atom in the fragment
      allocate(fc(nfrag,m3,m3))
c fc = the force constant matrix (hessian) of each atom in the fragment
      allocate(nat(nfrag,maxatom))
c nat = the whole-molecule atom number for real atoms
      allocate(numstore(nfrag,maxatom))
c numstore = the atomic number of each atom in each fragment
c eg 6 for carbon

c now rewind and read the data into the necessary arrays
      rewind(unit=5)
      read(5,*)
      read(5,*)
      do m=1,nfrag

      read(5,*)natomfrag(m)

      do n=1,natomfrag(m)
         read(5,*)ijk,jk,irub,(c(m,n,k),k=1,3)
           do i=1,3
           c(m,n,i)=c(m,n,i)/bohr
           enddo
       numstore(m,n)=jk
      enddo

      read(5,*)energyfrag(m)
c     read(5,*)comlin
c     read(5,*)comlin
c     read(5,*)comlin

      if(nflag.ge.1)then
      do n=1,natomfrag(m)
      read(5,*)ijk,ijk,(force(m,n,k),k=1,3)
      enddo
      else
      do n=1,natomfrag(m)
      do k=1,3
       force(m,n,k)=0.d0
      enddo
      enddo
      endif
c     read(5,*)comlin 

      if(nflag.eq.2)then

      n3=3*natomfrag(m)

      do i=1,n3
      do j=1,i
       read(5,*)fc(m,i,j)
      enddo
      enddo
c symmetrize fc
      do i=1,n3
      do j=i,n3
      fc(m,i,j)=fc(m,j,i)
      enddo
      enddo

      else
c nflag=1
      do i=1,n3
      do j=1,n3
      fc(m,i,j)=0.d0
      enddo
      enddo
c end the nflag if
      endif


c  end the m loop over frags
      enddo

      
      return
      end

