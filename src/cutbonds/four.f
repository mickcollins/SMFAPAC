      subroutine four

      use fractdata
      implicit double precision(a-h,o-z)

c intended to identify 4 member rings and to
c amend the ma matrix to ensure that such rings cannot be
c broken in such a way as to lead to interacting caps

c arrays ifour, ma, ilink are in common

c     dimension ifour(natomm,natomm)
      dimension ifour(500,4),ifr(20),ifr0(20)

      dimension mneg(12)

      ifour=0

      nafrag=numat(nfrag)
      if(nafrag.lt.4)return

      nring=0
c nearbours in a 4 ring?
      do n1=1,nafrag
       neg=0
       if(itp(n1).eq.-1)go to 91
      do k1=1,itp(n1)+1
       if(ibo(n1,k1).lt.0)neg=neg+1
       n2=iabs(ibo(n1,k1))
       if(itp(n2).eq.-1)go to 92
       do k2=1,itp(n2)+1
       do k3=1,itp(n1)+1
        ineg=0
        if(iabs(ibo(n2,k2)).eq.n1)go to 9
        if(iabs(ibo(n1,k3)).eq.n2)go to 9
        if(iabs(ilink(iabs(ibo(n2,k2)),iabs(ibo(n1,k3)))).eq.1)then
         if(ibo(n2,k2).lt.0)ineg=ineg+1
         if(ibo(n1,k3).lt.0)ineg=ineg+1
         if(neg+ineg.ge.3)go to 9
c we have a 4 ring
        ifr(1)=n1
        ifr(2)=n2
        ifr(3)=iabs(ibo(n2,k2))
        ifr(4)=iabs(ibo(n1,k3))
        call piksrt(4,ifr)
         if(nring.gt.0)then
          do k=1,nring
           idiff=0
           do j=1,4
           idiff=idiff+iabs(ifour(k,j)-ifr(j))
           enddo
           if(idiff.eq.0)go to 9
          enddo
         endif
        nring=nring+1
c       write(6,*)' ring in first part of four'
c       write(6,*)n1,ibo(n1,k1),ibo(n2,k2),ibo(n1,k3)
         do j=1,4
          ifour(nring,j)=ifr(j)
         enddo
c        ifour(n1,n2)=1
c        ifour(n2,n1)=1
        endif
9     continue
       enddo
       enddo
92      enddo
91      enddo

c opposite in a 4 ring?

      do n1=1,nafrag-1
      do n2=n1+1,nafrag
       match=0
       do i=1,12
        mneg(i)=0
       enddo
       if(itp(n1).eq.-1)go to 93
       do k1=1,itp(n1)+1
       if(itp(n2).eq.-1)go to 94
       do k2=1,itp(n2)+1
 
        if(iabs(ibo(n1,k1)).eq.iabs(ibo(n2,k2)))then
         match=match+1
         if(ibo(n1,k1).lt.0)mneg(match)=mneg(match)+1
         if(ibo(n2,k2).lt.0)mneg(match)=mneg(match)+1
         ifr(2+match)=iabs(ibo(n1,k1))
        endif
       enddo
       enddo
       if(match.ge.2)then

c we have one or more 4 ring
        do j0=1,match
         ifr0(j0)=ifr(2+j0)
        enddo
        do j1=1,match-1
        do j2=j1+1,match

          if(mneg(j1)+mneg(j2).ge.3)go to 99

        ifr(1)=n1
        ifr(2)=n2
        ifr(3)=ifr0(j1)
        ifr(4)=ifr0(j2)

        call piksrt(4,ifr)
         if(nring.gt.0)then
          do k=1,nring
           idiff=0
           do j=1,4
           idiff=idiff+iabs(ifour(k,j)-ifr(j))
           enddo
           if(idiff.eq.0)go to 99
          enddo
c end the ring if
         endif
          nring=nring+1
c       write(6,*)' ring in second part of four'
          do j=1,4
           ifour(nring,j)=ifr(j)
          enddo
99     continue
c end the j1,j2 loops
       enddo
       enddo
c end the match if
       endif

94      enddo
93      enddo

c     write(6,*)' four has nring = ',nring
c     write(6,*)'nfrag = ',nfrag,nafrag
c      if(nring.gt.0)then
c        write(6,*)' ring'
c        do k=1,nring
c        write(6,*)(ifour(k,j),j=1,4)
c        enddo
c      endif

c     do n1=1,nafrag
c     do n2=1,nafrag
c      if(ifour(n1,n2).eq.1)write(6,*)n1,n2
c     enddo
c     enddo

      if(nring.eq.0)return

c now amend the ma array

      jcount=0
100   ncount=0
      do n=1,nafrag
       do i=1,nring
        isum=0
        do j=1,4
        isum=isum+ma(n,ifour(i,j))
        enddo
        if(isum.gt.1.and.isum.lt.4)then
         ncount=ncount+1
         do j=1,4
          ma(n,ifour(i,j))=1
          ma(ifour(i,j),n)=1
         enddo
        endif
       enddo
      enddo

c debug in four.f
c     if(nafrag.eq.6)then
c      write(6,*)' ifour values'
c      do i=1,nring
c      write(6,*)(ifour(i,j),j=1,4)
c      enddo
c      write(6,*)' ma in four'
c      do n=1,nafrag
c       write(6,2222)(ma(n,k),k=1,nafrag)
c      enddo
c      write(6,*)
c2222   format(8i4)
c     stop
c     endif

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
c         if(ifour(n1,n2).eq.1.and.ma(n,n2).eq.1)match=1
c11      continue
c        if(match.eq.1)then
c         do n3=1,nafrag
c          if(ifour(n1,n3).eq.1)ma(n,n3)=1
c         enddo
c        endif
c10     continue
c      enddo

      return
      end


