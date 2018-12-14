      subroutine cancel

      use fractdata
      implicit double precision(a-h,o-z)

      dimension ibond1(nsmall,nsmall),ibond2(nsmall,nsmall)
c     dimension lki(natomm,natomm),lkj(natomm,natomm)

c get rid of cancelling fragments

      do 70 i=1,nf-1
      if(nstop(i).eq.2)go to 70

       ic=0
       do n=1,numat(i)
       do m=1,itype(i,n)+1
        ic=ic+1
        ibond1(n,m)=ib1(i,ic)
       enddo
       enddo

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

       ic=0
       do n=1,numat(j)
       do m=1,itype(j,n)+1
        ic=ic+1
        ibond2(n,m)=ib1(j,ic)
       enddo
       enddo

c compare atoms
      match=0
      do n1=1,numat(i)
      do n2=1,numat(j)
       if(natstore(i,n1).eq.natstore(j,n2))then
c compare bonding
        if(itype(i,n1).ne.itype(j,n2))go to 71
        like=0
        do k1=1,itype(i,n1)+1
        do k2=1,itype(j,n2)+1
         if(natstore(i,iabs(ibond1(n1,k1))).eq.natstore(j,iabs(ibond2(n2,k2))))then
         like=like+1
         go to 7100
         endif
        enddo
7100    enddo
        if(like.eq.itype(i,n1)+1)then
        match=match+1
        go to 710
        endif
       endif
      enddo
710   enddo
      if(match.eq.numat(i))then
       nstop(i)=2
       nstop(j)=2
       go to 70
      endif

71    continue
70    continue

c  reset the arrays

      ic=0
      do k=1,nf
      if(nstop(k).le.1)then
      ic=ic+1
       isign(ic)=isign(k)
      nstop(ic)=nstop(k)
      do i=1,nsmall
      itype(ic,i)=itype(k,i)
      enddo
c     do i=1,nsmall
c     do j=1,nsmall
c     ibond(ic,i,j)=ibond(k,i,j)
c     enddo
c     enddo

      do j=1,3*nsmall
       ib1(ic,j)=ib1(k,j)
      enddo

      numat(ic)=numat(k)
      do i=1,numat(k)
      natstore(ic,i)=natstore(k,i)
      enddo
      endif
      enddo

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
      do i=1,3*nsmall
       ib1(jc,i)=0
      enddo

      do i=1,nsmall
      natstore(jc,i)=0
      enddo
      numat(jc)=0

      enddo

      idiff=nf-ic
      nf=ic
      nfrag=1
c     write(6,*)nfrag,nf
      write(6,*)

      return
      end
