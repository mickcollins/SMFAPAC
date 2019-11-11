      subroutine expand(iter)

      use fractdata
      implicit double precision(a-h,o-z)


c iter is the expansion level
c we expand the icn type into the link, ibond etc

      dimension natsf(nsmall),nlink(nsmall,nsmall),nold(nsmall)

c  expand each fragment in turn


      do i=1,nf

       do n=1,nsmall
       do m=1,nsmall
        ibond(n,m)=0
       enddo
       enddo

       ic=0
       do n=1,numat(i)
       do m=1,itype(i,n)+1
       ic=ic+1
       ibond(n,m)=ib1(i,ic)
       enddo
       enddo

       n1=0
       do n=1,numat(i)
         n1=n1+1
         if(n1.gt.nsmall)then
          write(6,*)' The nsmall parameter is too small in expand'
          stop
         endif
         natsf(n1)=ngpf(natstore(i,n),1,iter)
         nold(n1)=n
         if(ngpf(natstore(i,n),2,iter).gt.0)then
         n1=n1+1
         if(n1.gt.nsmall)then
          write(6,*)' The nsmall parameter is too small in expand'
          stop
         endif
         natsf(n1)=ngpf(natstore(i,n),2,iter)
         nold(n1)=n
         endif
       enddo

       nnew=n1

c record a link array of the current bonding
       do j1=1,nsmall
       do j2=1,nsmall
        nlink(j1,j2)=0
       enddo
       enddo

      do j1=1,nnew
       n=nold(j1)
       do k1=1,itype(i,n)+1
        l=iabs(ibond(n,k1))
        do j2=1,nnew
         if(l.eq.nold(j2))then
          nlink(j1,j2)=1
          nlink(j2,j1)=1
         endif
        enddo
       enddo
      enddo

      do j1=1,nnew
      do j2=1,nnew
       if(nold(j1).eq.nold(j2))then
        nlink(j1,j2)=1
        nlink(j2,j1)=1
       endif
      enddo
      enddo

c nlink ne 0 only if j1 and j2 are part of groups
c bonding at the old iteration

c now re-constitute the fragment

      numat(i) = nnew

      do n=1,nsmall
       natstore(i,n)=0
       itype(i,n)=-1
       do m=1,nsmall
        ibond(n,m)=0
       enddo
      enddo

      do n=1,numat(i)
       natstore(i,n)=natsf(n)
      enddo


      do n=1,numat(i)
       ic=0
       do k=1,itf(natstore(i,n),iter-1)+1
c check that all neighbors are actually in this fragment
        do m=1,numat(i)
         mm=ibf(natstore(i,n),k,iter-1)
         if(iabs(mm).eq.natstore(i,m))then
          if(nlink(n,m).gt.0)then
          ic=ic+1
          ibond(n,ic)=m*(mm/iabs(mm))
          go to 100
         endif
         endif
        enddo
100    continue
       enddo
c       itype(i,natstore(i,n))=ic-1
        itype(i,n)=ic-1
c end loop over n
      enddo

       ic=0
       do n=1,numat(i)
       do m=1,itype(i,n)+1
       ic=ic+1
       if(ic.gt.6*nsmall)then
        write(6,*)' bonds exceed 6*nsmall in subroutine expand'
        stop
       endif
       ib1(i,ic)=ibond(n,m)
       enddo
       enddo

      it=0
      do n=1,numat(i)
      if(itype(i,n).gt.it)it=itype(i,n)
      enddo
      if(numat(i).le.it+Level)then
      nstop(i)=1
      else
      nstop(i)=0
      endif

c     do n=numat(i)+1,nsmall
c     do k=1,nsmall
c      ibond(i,n,k)=0
c     enddo
c     enddo
c     do n=1,numat(i)
c     do k=itype(i,n)+2,nsmall
c      ibond(i,n,k)=0
c     enddo
c     enddo

c end loop over fragments
      enddo

      natom=numgroups(iter-1)

      nfrag=1
      write(6,*)' After expansion there are ',natom,' groups'

      write(66,*)' After expansion'
      write(66,*)
      write(66,*)

      write(66,*)' The fragments and bonding are'
      do i=1,nf
       ic=0
       do n=1,numat(i)
       do m=1,itype(i,n)+1
       ic=ic+1
       ibond(n,m)=ib1(i,ic)
       enddo
       enddo
      write(66,*)' fragment ',i
      write(66,*)' atoms ',(natstore(i,k),k=1,numat(i))
      write(66,*)' connectivity '
      do m=1,numat(i)
       write(66,*)m,itype(i,m),(ibond(m,k),k=1,itype(i,m)+1)
      enddo
      enddo
      write(66,*)
      write(66,*)' The fragments in natural labels and bonding are'
      do i=1,nf
       ic=0
       do n=1,numat(i)
       do m=1,itype(i,n)+1
       ic=ic+1
       ibond(n,m)=ib1(i,ic)
       enddo
       enddo
      write(66,*)' fragment ',i
      write(66,*)' atoms ',(natstore(i,k),k=1,numat(i))
      write(66,*)' connectivity '
      do m=1,numat(i)
       write(66,*)natstore(i,m),itype(i,m),
     .          (natstore(i,iabs(ibond(m,k))),k=1,itype(i,m)+1)
      enddo
      enddo
      write(66,*)


      return
      end

