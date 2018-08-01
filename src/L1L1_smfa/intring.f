      subroutine intring(ngroups,nom,k1,k2,ik1k2
     .             ,ic,matstore,nsign,nffinal0,nffinal1,
     .             numat,natstore,nffinal,nbextra)

      implicit double precision(a-h,o-z)

      dimension nom(ngroups,ngroups),
c    .          matstore(10*nffinal0,2),nsign(10*nffinal0)
     .          matstore(nffinal0,2),nsign(nffinal0)


      dimension natstorek1(8),natstorek2(8),mata(8,8)
      dimension nk1(8),nk2(8)

      dimension numat(2*nffinal1),natstore(2*nffinal1,8)

c this subroutine evaluates the L1--L1 interaction
c where at least one fragment is a ring

c rings cannot be broken without leading to closely
c spaced caps. So if the L1..L1 interaction is fully
c allowed, it is implemented in full
c If totally disallowed, then no result
c If partially allowed, we break the ring into groups

c generate a local exclusion matrix

      numatk1=numat(k1)
      numatk2=numat(k2)
      do j=1,8
       natstorek1(j)=natstore(k1,j)
       natstorek2(j)=natstore(k2,j)
      enddo


c     write(6,*)numatk1,numatk2
c     stop
      ie=0
      id=0
      do n1=1,numatk1
      do n2=1,numatk2
       ie=ie+1
       mata(n1,n2)=nom(natstorek1(n1),natstorek2(n2))
       id=id+mata(n1,n2)
      enddo
      enddo
c     do n1=1,numatk1
c      write(6,*)(nom(natstorek1(n1),natstorek2(k)),k=1,numatk2)
c     enddo
c     stop

      if(id.eq.0)then
c nothing is allowed, so nothing to do
       return
      endif

      if(id.eq.ie)then
c all allowed
       if(ic.ge.10*nffinal0)then
        write(6,*)' ic too big'
      call cancelL1L1(10*nffinal0,ic,matstore,nsign)
       endif

       ic=ic+1
       matstore(ic,1)=k1
       matstore(ic,2)=k2
       nsign(ic)=ik1k2
       return
      endif

c see if any full fragment calcs are possible

      if(numatk2.gt.numatk1)then

       do n1=1,numatk1
        nk1(n1)=0
       do n2=1,numatk2
        nk1(n1)=nk1(n1)+mata(n1,n2)
       enddo
       enddo

       do n1=1,numatk1
        if(nk1(n1).eq.numatk2)then
       if(ic.ge.10*nffinal0)then
        write(6,*)' ic too big'
      call cancelL1L1(10*nffinal0,ic,matstore,nsign)
       endif

        ic=ic+1
        matstore(ic,2)=k2
        mg=natstorek1(n1)
      call match_mg(2*nffinal1,mg,numat,natstore,nffinal,nbextra,match)
c complete the interaction
        matstore(ic,1)=match
        nsign(ic)=ik1k2
c wipe out these mata entries
        do n2=1,numatk2
         mata(n1,n2)=0
        enddo
c end the if(nk1(n1).eq.numatk2)
       endif

c end loop over n1
       enddo

c else for       if(numatk2.gt.numatk1)
      else

       do n2=1,numatk2
        nk2(n2)=0
       do n1=1,numatk1
        nk2(n2)=nk2(n2)+mata(n1,n2)
       enddo
       enddo

       do n2=1,numatk2

        if(nk2(n2).eq.numatk1)then

       if(ic.ge.10*nffinal0)then
        write(6,*)' ic too big'
      call cancelL1L1(10*nffinal0,ic,matstore,nsign)
       endif

        ic=ic+1
        matstore(ic,1)=k1
        mg=natstorek2(n2)
      call match_mg(2*nffinal1,mg,numat,natstore,nffinal,nbextra,match)
c complete the interaction
        matstore(ic,2)=match
        nsign(ic)=ik1k2
c wipe out these mata entries
        do n1=1,numatk1
         mata(n1,n2)=0
        enddo
c end the if(nk2(n2).eq.numatk1)
       endif

c end loop over n1
       enddo

c end if(numatk2.gt.numatk1)
      endif

c any remaining mata not zero are done by group

      do n1=1,numatk1
      do n2=1,numatk2
       if(mata(n1,n2).ne.0)then

       mg=natstorek1(n1)
      call match_mg(2*nffinal1,mg,numat,natstore,nffinal,nbextra,match)
c save match
       match1=match

       mg=natstorek2(n2)
      call match_mg(2*nffinal1,mg,numat,natstore,nffinal,nbextra,match)
c save match
       match2=match
c record the fragment numbers and the net sign
       if(ic.ge.10*nffinal0)then
        write(6,*)' ic too big'
      call cancelL1L1(10*nffinal0,ic,matstore,nsign)
       endif

       ic=ic+1
       matstore(ic,1)=match1
       matstore(ic,2)=match2
       nsign(ic)=ik1k2
      endif
c end the n1,n2 loops
      enddo
      enddo

      return
      end
