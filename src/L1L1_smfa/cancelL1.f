      subroutine cancelL1(ndim,nffinal,numat,natstore,ksign,nat1,nch1,lab1,c1)
      implicit double precision(a-h,o-z)

      dimension numat(ndim),ksign(ndim),natstore(ndim,8),
     .          nch1(ndim),nat1(ndim)
      dimension c1(ndim,60,3)

      character*2 lab1(ndim,60)

      do n1=1,nffinal-1
      if(ksign(n1).eq.0)go to 1
      do n2=n1+1,nffinal
      if(ksign(n2).eq.0)go to 2
      if(numat(n2).ne.numat(n1))go to 2
      do j=1,numat(n1)
       if(natstore(n1,j).ne.natstore(n2,j))go to 2
      enddo
c we have a match
      ksign(n1)=ksign(n1)+ksign(n2)
      ksign(n2)=0
2     continue
      enddo
1     continue
      enddo

c elimiate zeros

      ic=0
      do n=1,nffinal
       if(ksign(n).ne.0)then

        ic=ic+1
        numat(ic)=numat(n)
        do j=1,numat(n)
         natstore(ic,j)=natstore(n,j)
        enddo
        do j=1,nat1(n)
         lab1(ic,j)=lab1(n,j)
         do k=1,3
          c1(ic,j,k)=c1(n,j,k)
         enddo
        enddo
        nat1(ic)=nat1(n)
        ksign(ic)=ksign(n)
        nch1(ic)=nch1(n)
        
       endif
      enddo

      if(ic.lt.nffinal)then
       do n=ic+1,nffinal
        numat(n)=0
        do j=1,8
         natstore(n,j)=0
        enddo
        do j=1,60
        do k=1,3
         c1(n,j,k)=0.d0
        enddo
        enddo
        nat1(n)=0
        ksign(n)=0
        nch1(n)=0
       enddo
      endif

      nffinal=ic

      return
      end
