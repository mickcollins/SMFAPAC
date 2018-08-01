      subroutine five

      use fractdata
      implicit double precision(a-h,o-z)

c intended to identify 5 member rings and to
c amend the ma matrix to ensure that such rings cannot be
c broken in such a way as to lead to interacting caps

c arrays ifive, ma, ilink are in common

c     dimension ifive(natomm,natomm)
      dimension ifive(500,5),ifr(5)

      ifive=0

      nafrag=numat(nfrag)
      if(nafrag.lt.5)return

      nring=0
c nearbours in a 5 ring?
      do n1=1,nafrag
      do k1=1,itp(n1)+1
       n2=ibo(n1,k1)
       do k2=1,itp(n2)+1
       do k3=1,itp(n1)+1
        if(ibo(n2,k2).eq.n1)go to 9
        if(ibo(n1,k3).eq.n2)go to 9
        n3=ibo(n2,k2)
        n4=ibo(n1,k3)
        if(n3.eq.n4)go to 9
         do k4=1,itp(n3)+1
         do k5=1,itp(n4)+1
        if(ibo(n3,k4).eq.n2)go to 99
        if(ibo(n3,k4).eq.n4)go to 99
        if(ibo(n4,k5).eq.n1)go to 99
        if(ibo(n4,k5).eq.n3)go to 99
         if(ibo(n3,k4).eq.ibo(n4,k5))then
c we have a 5 ring
           ifr(1)=n1
           ifr(2)=n2
           ifr(3)=n3
           ifr(4)=n4
           ifr(5)=ibo(n3,k4)
           call piksrt(5,ifr)
           if(nring.gt.0)then
            do k=1,nring
             idiff=0
             do j=1,5
              idiff=idiff+iabs(ifive(k,j)-ifr(j))
             enddo
            if(idiff.eq.0)go to 99
            enddo
c end the nring gt 0 if
           endif
           nring=nring+1
          do j=1,5
           ifive(nring,j)=ifr(j)
          enddo
c end the ibo match if
         endif
99    continue
        enddo
        enddo
9     continue
       enddo
       enddo
      enddo
      enddo

c     write(6,*)'nfrag = ',nfrag,nafrag
c      if(nring.gt.0)then
c        write(6,*)' ring'
c        do k=1,nring
c        write(6,*)(ifive(k,j),j=1,5)
c        enddo
c      endif

c     do n1=1,nafrag
c     do n2=1,nafrag
c      if(ifive(n1,n2).eq.1)write(6,*)n1,n2
c     enddo
c     enddo

      if(nring.eq.0)return

c now amend the ma array

100   ncount=0
      do n=1,nafrag
       do i=1,nring
        isum=0
        do j=1,5
        isum=isum+ma(n,ifive(i,j))
        enddo
        if(isum.gt.2.and.isum.lt.5)then
         ncount=ncount+1
         do j=1,5
          ma(n,ifive(i,j))=1
          ma(ifive(i,j),n)=1
         enddo
        endif
       enddo
      enddo

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
c         if(ifive(n1,n2).eq.1.and.ma(n,n2).eq.1)match=1
c11      continue
c        if(match.eq.1)then
c         do n3=1,nafrag
c          if(ifive(n1,n3).eq.1)ma(n,n3)=1
c         enddo
c        endif
c10     continue
c      enddo

      return
      end


