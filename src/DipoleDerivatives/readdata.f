      subroutine readdata
      use derivheader

c the force and (possibly) hessian data for each fragment is 
c read from standard input 


      bohr=1.d0/1.8897259886d0

      open(unit=1,file='IN_PACKAGE',status='old')
       read(1,*)npack
      close(unit=1)

c now read the data
      read(5,*)nfrag
c nfrag = the number of fragments

c fact = the factor which allocates the fraction of a cap force
c to each of the two real atoms

c now have the max number of atoms in a fragment
      m3=3*maxatom

      allocate(dipdr(3,3,natom))
      allocate(dipdrfrag(3,3,maxatom))
      write(6,*)
      write(6,*)nfrag
      write(6,*)
      do m=1,nfrag

c atoms in the frag
       read(5,*)numatom
       if(numatom.ne.natomfrag(m))then
       write(6,*)'Number of atoms in FragDipderivs,',numatom,' does not'
        write(6,*)'equal the number of atom in frags.out for'
        write(6,*)'fragment ',m
        stop
       endif
c first read the data for each fragment, according to the
c format of the particular qc program
c the first index refers to the xyz component of the dipole
c the second index refers to the xyz of the atom
c the third index refers to the atom number
       if(npack.eq.1)then
        do n=1,numatom
         do k=1,3
          read(5,*)(dipdrfrag(j,k,n),j=1,3)
         enddo
        enddo
       endif
       if(npack.eq.2)then
        read(5,*)(((dipdrfrag(j,k,n),j=1,3),k=1,3),n=1,numatom)
       endif
       if(npack.eq.3)then
        do j=1,3
         do n=1,numatom
          do k=1,3
           read(5,*)dipdrfrag(j,k,n)
          enddo
         enddo
        enddo
       endif
       if(npack.eq.4)then
        do n=1,numatom
         do k=1,3
          read(5,*)(dipdrfrag(j,k,n),j=1,3)
         enddo
        enddo
       endif

c now assign to the real atoms
       do n=1,nat0(m)
       do k=1,3
        do j=1,3
         dipdr(j,k,nat(m,n))=dipdr(j,k,nat(m,n))
     .                       +dipdrfrag(j,k,n)*isign(m)
        enddo
       enddo
       enddo

       do n=1,ncap(m)
       do k=1,3
        do j=1,3
        dipdr(j,k,natcap(m,n,1))=dipdr(j,k,natcap(m,n,1))+
     .              (1.d0-fact(m,n))*dipdrfrag(j,k,nat0(m)+n)*isign(m)
        dipdr(j,k,natcap(m,n,2))=dipdr(j,k,natcap(m,n,2))+
     .               fact(m,n)*dipdrfrag(j,k,nat0(m)+n)*isign(m)
        enddo
       enddo
       enddo
c end loop over frags
      enddo
      
      return
      end

