      program combfreq

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

c readin the forces and hessians of the fragments
c some necessary arrays are also defined 
      call readdata

c read in the atom numbering of the fragments
      call readfrag

c assign the mass, elemental symbol and cartsian coords
c for each atom in the whole molecule
      call massset


c read in the coefficients for each fragment in the energy sum
      call readsigns

c zero arrays

      energy=0.d0

      allocate(dvdx(natom,3))
      allocate(fa(3*natom,3*natom))
c dvdx = the gradient for the whole molecule
c fa = the hessian for the whole molecule

      do n1=1,natom
      do k1=1,3
      i1=3*(n1-1)+k1
      dvdx(n1,k1)=0.d0
      do n2=1,natom
      do k2=1,3
      i2=3*(n2-1)+k2
      fa(i1,i2)=0.d0
      enddo
      enddo
      enddo
      enddo

c  combine the energies

c     write(6,*)energy
      do m =1,nfrag
      write(6,*)energyfrag(m)
      energy=energy+dble(float(isign(m)))*energyfrag(m)
      enddo
c     write(6,33)energy
      
c  combine the first derivatives

      do m =1,nfrag

      do n=1,nat0(m)
      do k=1,3
      dvdx(nat(m,n),k)=dvdx(nat(m,n),k)-isign(m)*force(m,n,k)
c the - sign in "-isign" converts forces to gradients
      enddo
      enddo
      if(ncap(m).eq.0)go to 111
c allocate forces on caps to the associated real atoms
      do n=1,ncap(m)
      do k=1,3
      dvdx(natcap(m,n,1),k)=dvdx(natcap(m,n,1),k)-isign(m)*
     .   (1.d0-fact(m,n))*force(m,nat0(m)+n,k)
      dvdx(natcap(m,n,2),k)=dvdx(natcap(m,n,2),k)-isign(m)*
     .   fact(m,n)*force(m,nat0(m)+n,k)
      enddo
      enddo
111   continue

      enddo


      if(nflag.eq.1)go to 222

c  combine the second derivatives

c  check symmetry
      do m =1,nfrag
      do n1=1,nat0(m)
      do k1=1,3
      i1=3*(n1-1)+k1
      do n2=1,nat0(m)
      do k2=1,3
      i2=3*(n2-1)+k2
      diff=fc(m,i1,i2)-fc(m,i2,i1)
      if(abs(diff).gt.1.d-12)then
      write(6,*)' fragment hessian not symmetric in Derivs'
      write(6,*)m,i1,i2,diff
      stop
      endif
      enddo
      enddo
      enddo
      enddo

      enddo

c allocate the second derivatives from the "real" atoms
      do m =1,nfrag

      do n1=1,nat0(m)
      do k1=1,3
      i1=3*(n1-1)+k1
      j1=3*(nat(m,n1)-1)+k1

      do n2=1,nat0(m)
      do k2=1,3
      i2=3*(n2-1)+k2
      j2=3*(nat(m,n2)-1)+k2

      fa(j1,j2)=fa(j1,j2)+isign(m)*fc(m,i1,i2)

      enddo
      enddo

      enddo
      enddo

c allocate the cross derivatives for a real atom and a cap
      do n1=1,nat0(m)
      do k1=1,3
      i1=3*(n1-1)+k1
      j1=3*(nat(m,n1)-1)+k1

      if(ncap(m).eq.0)go to 112
      do n2=1,ncap(m)
      do k2=1,3
      i2=3*(nat0(m)+n2-1)+k2
      ja2=3*(natcap(m,n2,1)-1)+k2
      jb2=3*(natcap(m,n2,2)-1)+k2
      fa(j1,ja2)=fa(j1,ja2)+isign(m)*(1.d0-fact(m,n2))*fc(m,i1,i2)
      fa(j1,jb2)=fa(j1,jb2)+isign(m)*fact(m,n2)*fc(m,i1,i2)
      enddo
      enddo
112   continue
      enddo
      enddo

c allocate the cross derivatives for a cap and a real atom
      if(ncap(m).eq.0)go to 113
      do n2=1,ncap(m)
      do k2=1,3
      i2=3*(nat0(m)+n2-1)+k2
      ja2=3*(natcap(m,n2,1)-1)+k2
      jb2=3*(natcap(m,n2,2)-1)+k2

      do n1=1,nat0(m)
      do k1=1,3
      i1=3*(n1-1)+k1
      j1=3*(nat(m,n1)-1)+k1

      fa(ja2,j1)=fa(ja2,j1)+isign(m)*(1.d0-fact(m,n2))*fc(m,i2,i1)
      fa(jb2,j1)=fa(jb2,j1)+isign(m)*fact(m,n2)*fc(m,i2,i1)

      enddo
      enddo

      enddo
      enddo

c allocate the derivatives from two caps
      do n2=1,ncap(m)
      do k2=1,3
      i2=3*(nat0(m)+n2-1)+k2
      ja2=3*(natcap(m,n2,1)-1)+k2
      jb2=3*(natcap(m,n2,2)-1)+k2

      do n1=1,ncap(m)
      do k1=1,3
      i1=3*(nat0(m)+n1-1)+k1
      ja1=3*(natcap(m,n1,1)-1)+k1
      jb1=3*(natcap(m,n1,2)-1)+k1

      fa(ja2,ja1)=fa(ja2,ja1)+isign(m)*
     . (1.d0-fact(m,n2))*(1.d0-fact(m,n1))*fc(m,i2,i1)
      fa(ja2,jb1)=fa(ja2,jb1)+isign(m)*
     . (1.d0-fact(m,n2))*fact(m,n1)*fc(m,i2,i1)
      fa(jb2,ja1)=fa(jb2,ja1)+isign(m)*
     . fact(m,n2)*(1.d0-fact(m,n1))*fc(m,i2,i1)
      fa(jb2,jb1)=fa(jb2,jb1)+isign(m)*
     . fact(m,n2)*fact(m,n1)*fc(m,i2,i1)

      enddo
      enddo

      enddo
      enddo
113   continue

c  end the mloop over nfrag
      enddo

c  check symmetyry
      do i1=1,3*natom
      do i2=i1,3*natom
      diff=fa(i1,i2)-fa(i2,i1)
      if(abs(diff).gt.1.d-12)then
       write(6,*)' Total hessian is not symmetric in Derivs'
       stop
      endif
      enddo
      enddo

222   continue

c  write out the derivatives, with masses, elemental symbols
c and coordinates
      open(unit=1,file='combinedFRAGderivs',status='unknown')

      write(1,*)" The flag is"
      write(1,*)nflag
      write(1,*)" The number of atoms is"
      write(1,*)natom
      write(1,*)" The masses are "
      do n=1,natom
      write(1,*)amas(n)
      enddo
      write(1,*)" The elements are "
      do n=1,natom
      write(1,50)lab(n)
      enddo
50    format(a2)
      write(1,*)" The energy is"
      write(1,33)energy
33      format(1x,3D25.16)
      write(1,*)" The coordinates are"
      do n=1,natom
      write(1,33)(coord(n,k),k=1,3)
      enddo
      write(1,*)" First derivatives"
      do n=1,natom
      write(1,33)dvdx(n,1)
      write(1,33)dvdx(n,2)
      write(1,33)dvdx(n,3)
      write(33,*)(dvdx(n,k),k=1,3)
      enddo
      write(1,*)" Upper triangle of the second derivatives"
      n3=3*natom
100   format(6f10.6)
        write(1,100)((fa(i,j),j=1,i),i=1,n3)
      close(unit=1)

      end
