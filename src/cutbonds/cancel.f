      subroutine cancel

      use fractdata
      implicit double precision(a-h,o-z)

c     dimension lki(natomm,natomm),lkj(natomm,natomm)

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
      do n2=1,numat(j)
       if(natstore(i,n1).eq.natstore(j,n2))then
c compare bonding
        if(itype(i,n1).ne.itype(j,n2))go to 71
        like=0
        do k1=1,itype(i,n1)+1
        do k2=1,itype(j,n2)+1
         if(natstore(i,iabs(ibond(i,n1,k1))).eq.natstore(j,iabs(ibond(j,n2,k2))))then
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
      do i=1,nsmall
      do j=1,nsmall
      ibond(ic,i,j)=ibond(k,i,j)
      enddo
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
      do i=1,nsmall
      do j=1,nsmall
      ibond(jc,i,j)=0
      enddo
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
