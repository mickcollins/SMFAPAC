      subroutine finalcancel

      use fractdata
      implicit double precision(a-h,o-z)

      integer, allocatable  :: jnstop(:),jisign(:), jnumat(:)
      integer, allocatable  :: jitype(:,:), jib1(:,:),jnatstore(:,:)

c get rid of cancelling fragments

      do 70 i=1,nf-1
      if(nstop(i).eq.2)go to 70
      do 71 j=i+1,nf
      if(nstop(j).eq.2)go to 71
      if(iabs(isign(i)*numat(i)+isign(j)*numat(j)).gt.0)go to 71
c patch for 1 member fragments
      if(numat(i).eq.1.and.iabs(natstore(i,1)-natstore(j,1)).gt.0)then
       go to 71
      endif

      if(numat(i).eq.1.and.iabs(natstore(i,1)-natstore(j,1)).eq.0)then
       nstop(i)=2
       nstop(j)=2
       go to 70
      endif 

c compare atoms
      match=0
      do n1=1,numat(i)
      ig=1
      do n2=1,numat(j)
       if(natstore(i,n1).eq.natstore(j,n2))ig=0
      enddo
      match=match+ig
      enddo
      if(match.eq.0)then
       nstop(i)=2
       nstop(j)=2
       go to 70
      endif

71    continue
70    continue

      do i=1,nf
       if(numat(i).eq.0)nstop(i)=2
      enddo

c  reset the arrays
c Note: we put the +ve fragments first

      ic=0
      do k=1,nf
       if(nstop(k).le.1)then
        ic=ic+1
        isign(ic)=isign(k)
        nstop(ic)=nstop(k)
        do i=1,nsmall
         itype(ic,i)=itype(k,i)
        enddo
c       do i=1,nsmall
c       do j=1,nsmall
c        ibond(ic,i,j)=ibond(k,i,j)
c       enddo
c       enddo

        do j=1,6*nsmall
         ib1(ic,j)=ib1(k,j)
        enddo

        numat(ic)=numat(k)
        do i=1,numat(k)
         natstore(ic,i)=natstore(k,i)
        enddo
       endif
72    enddo


      do jc=ic+1,nfragm
      isign(jc)=1
      nstop(jc)=0
      do i=1,nsmall
      itype(jc,i)=-1
      enddo
c     do i=1,nsmall
c     do j=1,nsmall
c     ibond(jc,i,j)=0
c     enddo
c     enddo

      do j=1,6*nsmall
       ib1(jc,j)=0
      enddo

      do i=1,nsmall
      natstore(jc,i)=0
      enddo
      numat(jc)=0

      enddo

      idiff=nf-ic
      nf=ic
      nfrag=1
      write(6,*)

c put the positive fragments first

      allocate(jnstop(nf))
      allocate(jisign(nf))
      allocate(jnumat(nf))
      allocate(jitype(nf,nsmall))
      allocate(jib1(nf,6*nsmall))
      allocate(jnatstore(nf,nsmall))

      ic=0
      do k=1,nf
      if(isign(k).gt.0)then
       ic=ic+1
       jnstop(ic)=nstop(k)
       jisign(ic)=isign(k)
       do i=1,nsmall
        jitype(ic,i)=itype(k,i)
       enddo
       do i=1,6*nsmall
        jib1(ic,i)=ib1(k,i)
       enddo
       jnumat(ic)=numat(k)
       do i=1,numat(k)
        jnatstore(ic,i)=natstore(k,i)
       enddo
      endif
      enddo
      do k=1,nf
      if(isign(k).lt.0)then
       ic=ic+1
       jnstop(ic)=nstop(k)
       jisign(ic)=isign(k)
       do i=1,nsmall
        jitype(ic,i)=itype(k,i)
       enddo
       do i=1,6*nsmall
        jib1(ic,i)=ib1(k,i)
       enddo
       jnumat(ic)=numat(k)
       do i=1,numat(k)
        jnatstore(ic,i)=natstore(k,i)
       enddo
      endif
      enddo

      do k=1,nf

       isign(k)=jisign(k)
       nstop(k)=jnstop(k)
       do i=1,nsmall
        itype(k,i)=jitype(k,i)
       enddo
       do i=1,6*nsmall
        ib1(k,i)=jib1(k,i)
       enddo
       numat(k)=jnumat(k)
       do i=1,numat(k)
        natstore(k,i)=jnatstore(k,i)
       enddo

      enddo

      deallocate(jnstop)
      deallocate(jisign)
      deallocate(jnumat)
      deallocate(jitype)
      deallocate(jib1)
      deallocate(jnatstore)


      return
      end
