      subroutine writerawdata
      use fractdata
      implicit double precision(a-h,o-z)

c this subroutine is called in bigfract after fragsteps
c it writes the data about atoms, bonds and fragmnents
c that is needed to restart fragmentation

#ifdef __GFORTRAN__
      open(unit=1,file='RAWDATA',status='unknown')
#else
      open(unit=1,file='RAWDATA',status='unknown',buffered='YES')
#endif

      write(1,*)natomall
      do n=1,natomall
       write(1,*)atoms(n),ichg(n)
      enddo
      do n=1,natomall
       write(1,*)(c(n,k),k=1,3)
      enddo

      write(1,*)natom
      do n=1,natom
       write(1,*)nfam(n)
       write(1,*)(ifam(n,j),j=1,nfam(n))
      enddo

      write(1,*)nbondso
      do i=1,nbondso
       write(1,*)mbstore(i),nbstore(i),multstore(i)
      enddo

c now the fragments
      write(1,*)nf
      do k=1,nf
       write(1,*)nstop(k),isign(k),numat(k)
       write(1,*)(itype(k,i),i=1,numat(k))
       write(1,*)(natstore(k,i),i=1,numat(k))
       ic=0
       do i=1,numat(k)
       do j=1,itype(k,i)+1
       ic=ic+1
       ibond(i,j)=ib1(k,ic)
       enddo
       enddo
       do i=1,numat(k)
        if(itype(k,i).gt.-1)then
        write(1,*)(ibond(i,j),j=1,itype(k,i)+1)
        endif
       enddo
      enddo

      return
      end


