      subroutine setupfrags1

      use fractdata
      implicit double precision(a-h,o-z)


c after compress is finished, we copy the last
c iteration output into the ibond and itype arrays

c iterfinal is the last compression level reached

      natom=numgroups(iterfinal)

      write(6,*)' in setupfrags, natom = ',natom

      do n=1,nsmall
       natstore(1,n)=0
       itype(1,n)=-1
      do m=1,nsmall
c      ibond(1,n,m)=0
       ibond(n,m)=0
      enddo
      enddo

      n1=0
      do n=1,natom
       if(itf(n,iterfinal).eq.-1)go to 100
       n1=n1+1
       natstore(1,n1)=n
       itype(1,n1)=itf(n,iterfinal)
100   enddo

      numat(1)=n1

      do n=1,numat(1)
       do k=1,itype(1,n)+1
        do m=1,numat(1)
         ii=ibf(natstore(1,n),k,iterfinal)
         if(iabs(ibf(natstore(1,n),k,iterfinal)).eq.natstore(1,m))then
c         ibond(1,n,k)=m*(ii/iabs(ii))
          ibond(n,k)=m*(ii/iabs(ii))
         endif
        enddo
       enddo
      enddo

      i=0
      do n=1,numat(1)
       do k=1,itype(1,n)+1
        i=i+1
        ib1(1,i)=ibond(n,k)
       enddo
      enddo
       

      nfrag=1
      nf=1
      nstop(1)=0
      isign(1)=1

c account for isolated groups
      do n=1,natom
       if(itf(n,iterfinal).eq.-1)then
         nf=nf+1
         natstore(nf,1)=n
         itype(nf,1)=-1
c         ibond(nf,1,1)=0
          ibond(1,1)=0
          ib1(nf,1)=0
         numat(nf)=1
         isign(nf)=1
         nstop(nf)=1
       endif
      enddo
       

      do i=1,nf
      write(6,*)' An initial compressed fragment is'
      write(6,*)(natstore(i,n),n=1,numat(i))
      write(6,*)' The connectivity is'
      ic=0
      do n=1,numat(i)
c       write(6,*)n,itype(i,n),(ibond(i,n,k),k=1,itype(i,n)+1)
       write(6,*)n,itype(i,n),(ib1(i,k),k=ic+1,ic+itype(i,n)+1)
       ic=ic+itype(i,n)+1
      enddo
      enddo
      write(6,*)

 

      return
      end

