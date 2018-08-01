      subroutine three

      use fractdata
      implicit double precision(a-h,o-z)

c intended to identify 3 member rings and to
c amend the ma matrix to ensure that such rings cannot be
c broken in such a way as to lead to interacting caps

c arrays ithree, ma, ilink are in common

c     dimension ithree(natomm,natomm)
      dimension ithree(500,3),ifr(3)

      ithree=0

      nafrag=numat(nfrag)

      nring=0
c nearbours in a 3 ring?
      do n1=1,nafrag
      neg=0
      do k1=1,itp(n1)+1
       if(ibo(n1,k1).lt.0)neg=neg+1
       n2=iabs(ibo(n1,k1))
       do k2=1,itp(n2)+1
       do k3=1,itp(n1)+1
        ineg=0
        if(ibo(n2,k2).eq.n1)go to 9
        if(ibo(n1,k3).eq.n2)go to 9
        if(ibo(n2,k2).eq.ibo(n1,k3))then
         if(ibo(n2,k2).lt.0)ineg=ineg+1
         if(ibo(n1,k3).lt.0)ineg=ineg+1
         if(neg+ineg.ge.2)go to 9
c we have a 3 ring
        ifr(1)=n1
        ifr(2)=n2
        ifr(3)=ibo(n2,k2)
        call piksrt(3,ifr)
         if(nring.gt.0)then
          do k=1,nring
           idiff=0
           do j=1,3
           idiff=idiff+iabs(ithree(k,j)-ifr(j))
           enddo
           if(idiff.eq.0)go to 9
          enddo
         endif
        nring=nring+1
         do j=1,3
          ithree(nring,j)=ifr(j)
         enddo
        endif
9     continue
       enddo
       enddo
      enddo
      enddo

c      write(6,*)' three has nring = ',nring

      if(nring.eq.0)return

c now amend the ma array

      jcount=0
100   ncount=0
      do n=1,nafrag
       do i=1,nring
        isum=0
        do j=1,3
        isum=isum+ma(n,ithree(i,j))
        enddo
        if(isum.ge.1.and.isum.lt.3)then
         ncount=ncount+1
         do j=1,3
          ma(n,ithree(i,j))=1
          ma(ithree(i,j),n)=1
         enddo
        endif
       enddo
      enddo

      jcount=jcount+1
c     if(jcount.eq.2)return
c check if more needs to be done
      if(ncount.gt.0)go to 100

c      do n=1,nafrag
c
c       do 10 n1=1,nafrag
c        if(ma(n,n1).eq.0)go to 10
c
c        match=0
c        do 11 n2=1,nafrag
c         if(n2.eq.n1)go to 11
c         if(ithree(n1,n2).eq.1.and.ma(n,n2).eq.1)match=1
c11      continue
c        if(match.eq.1)then
c         do n3=1,nafrag
c          if(ithree(n1,n3).eq.1)ma(n,n3)=1
c         enddo
c        endif
c10     continue
c      enddo

      return
      end


