      subroutine compress

      use fractdata
      implicit double precision(a-h,o-z)

      integer, allocatable  :: ifin(:)
      allocate(ifin(natom))


c compress the initial set of groups

c first copy the original ibond and itype arrays


      numgroups(1) = natom
c     do n=1,natom
c     itf(n,1) = itype(1,n)
c     do k=1,itype(1,n)+1
c     ibf(n,k,1) = ibond(1,n,k)
c     enddo
c     enddo

c  now loop over the iterations

      do i=2,ncomp
c  find the pairs
      ng=0
      do n=1,numgroups(i-1)
       ifin(n)=0
      enddo
c first the isolated atoms
      do n=1,numgroups(i-1)
       if(itf(n,i-1).eq.-1)then
        write(6,*)'found -1 at ',i,n
        ng=ng+1
        ngpf(ng,1,i)=n
        ngpf(ng,2,i)=0
        ifin(n)=1
       endif
      enddo
c next the dangling atoms

      do n=1,numgroups(i-1)
      if(ifin(n).gt.0)go to 100
      if(itf(n,i-1).eq.0)then
       ifin(n)=1
       ng=ng+1
       ngpf(ng,1,i)=n
       if(ifin(iabs(ibf(n,1,i-1))).eq.0)then
        ngpf(ng,2,i)=iabs(ibf(n,1,i-1))
        ifin(iabs(ibf(n,1,i-1)))=1
       else
        ngpf(ng,2,i)=0
       endif
      endif
100   enddo

c  next the rest of the atoms
      do n=1,numgroups(i-1)
       if(ifin(n).gt.0)go to 101
c start a new group
       ng=ng+1
       ngpf(ng,1,i)=n
       ifin(n)=1
c try to add a second atom to the group from the neighbors of the first
       do k=1,itf(n,i-1)+1
        m=iabs(ibf(n,k,i-1))
        if(ifin(m).eq.1)go to 103
c check to see that this neighbour and n are not
c connected to a common atom
        match=0 
        do k1=1,itf(m,i-1)+1
        do k2=1,itf(n,i-1)+1
         if(iabs(ibf(n,k1,i-1)).eq.iabs(ibf(m,k2,i-1)))match=1
        enddo
        enddo
       if(match.eq.0)then
c n and m are connected to a common atom. Need to
c check these two atoms are not both connected to a next level group
        like=0
        if(ng.eq.1)go to 105
        do l=1,ng-1
        do k1=1,itf(m,i-1)+1
        do k2=1,itf(n,i-1)+1
         if(iabs(ibf(n,k1,i-1)).eq.ngpf(l,1,i).and.
     .      iabs(ibf(m,k2,i-1)).eq.ngpf(l,1,i))like=1
         if(iabs(ibf(n,k1,i-1)).eq.ngpf(l,2,i).and.
     .      iabs(ibf(m,k2,i-1)).eq.ngpf(l,2,i))like=1
         if(iabs(ibf(n,k1,i-1)).eq.ngpf(l,1,i).and.
     .      iabs(ibf(m,k2,i-1)).eq.ngpf(l,2,i))like=1
         if(iabs(ibf(n,k1,i-1)).eq.ngpf(l,2,i).and.
     .      iabs(ibf(m,k2,i-1)).eq.ngpf(l,1,i))like=1
         enddo
         enddo
         enddo

105     if(like.eq.0)then 
         ngpf(ng,2,i)=iabs(ibf(n,k,i-1))
         ifin(iabs(ibf(n,k,i-1)))=1
         go to 102
        endif
       endif

103    enddo
       ngpf(ng,2,i) = 0
102   continue
101   enddo

c  record how many pairs (new groups)
      numgroups(i) = ng
      write(6,*)' for compression i = ',i
      write(6,*)' numgroups(i) = ',numgroups(i)

c work out the new "ibond" and "itype"arrays

      do n1=1,ng
      j=0

      if(itf(ngpf(n1,1,i),i-1).eq.-1)then
       itf(n1,i)=-1
       ibf(n1,1,i)=0
       go to 304
      endif  

      do k1=1,itf(ngpf(n1,1,i),i-1)+1
       ja=ibf(ngpf(n1,1,i),k1,i-1)
       if(iabs(ja).eq.ngpf(n1,2,i))go to 201
        do n2=1,ng
         if(n2.eq.n1)go to 202
         if(iabs(ja).eq.ngpf(n2,1,i).or.iabs(ja).eq.ngpf(n2,2,i))then
          if(j.ge.1)then
           do j1=1,j
            if(iabs(ibf(n1,j1,i)).eq.n2)go to 202
           enddo
          endif
          j=j+1
          ibf(n1,j,i)=n2*(ja/iabs(ja))
         endif
202     enddo
201   enddo

      if(ngpf(n1,2,i).eq.0)go to 303

c repeat for ngpf(n1,2,i)
      do k1=1,itf(ngpf(n1,2,i),i-1)+1
       ja=ibf(ngpf(n1,2,i),k1,i-1)
       if(iabs(ja).eq.ngpf(n1,1,i))go to 301
        do n2=1,ng
         if(n2.eq.n1)go to 302
         
         if(iabs(ja).eq.ngpf(n2,1,i).or.iabs(ja).eq.ngpf(n2,2,i))then
c avoid double counting arising from triples and tetramers
         if(j.ge.1)then
         do j1=1,j
          if(iabs(ibf(n1,j1,i)).eq.n2)go to 302
         enddo
         endif
          j=j+1
          ibf(n1,j,i)=n2*(ja/iabs(ja))
         endif
302     enddo
301   enddo
303   continue
      itf(n1,i)=j-1
304   enddo


      write(6,*)' compression at'
       write(6,*)i,numgroups(i)
      write(6,*)' pairs'
      do n=1,numgroups(i)
      write(6,*)n,ngpf(n,1,i),ngpf(n,2,i)
      enddo
      write(6,*)' n, itf, ibf'
      do n=1,numgroups(i)
      write(6,*)n,itf(n,i),(ibf(n,k,i),k=1,itf(n,i)+1)
      enddo

c  finished with iteration i

c if all the new "pairs" are singles, we can stop at
c the previous iteration
      isum=0
      do j=1,numgroups(i)
      isum=isum+iabs(ngpf(j,2,i))
      enddo
      if(isum.eq.0)then
       iterfinal=i-1
       write(6,*)' Stopped at compression ',iterfinal
       write(6,*)' as all new pairs are singles'
       go to 104
      endif

c if the number of groups is below 10 we can stop here
       if(ng.le.10)then
       iterfinal=i
       write(6,*)' Stopped at compression ',iterfinal
       write(6,*)' as the number of groups is below 10'
       go to 104
       endif

c if the number of groups is not reducing, then we stop
c      if(numgroups(i).gt.numgroups(i-1)-5)then
c      iterfinal=i
c      write(6,*)' Stopped at compression ',iterfinal
c      write(6,*)' as the number of groups is not reducing'
c      go to 104
c      endif

c if the number of groups is not reducing, then we stop
       s1=dble(float(numgroups(i)))
       s2=dble(float(numgroups(i-1)))
       if((s1/s2).gt.0.98d0)then
       iterfinal=i
       write(6,*)' Stopped at compression ',iterfinal
       write(6,*)' as the number of groups is not reducing'
       go to 104
       endif

      if(i.eq.ncomp)then
      iterfinal=i
       write(6,*)' Stopped at compression ',iterfinal
       write(6,*)' as this is the limit given by parameter ncomp'
      go to 104
      endif

      enddo


104   continue


       write(6,*)' bottom of compress'
      write(6,*)i,numgroups(i)
      write(6,*)' pairs'
      do n=1,numgroups(i)
      write(6,*)n,ngpf(n,1,i),ngpf(n,2,i)
      enddo
      write(6,*)' n, itf, ibf'
      do n=1,numgroups(i)
      write(6,*)n,itf(n,i),(ibf(n,k,i),k=1,itf(n,i)+1)
      enddo

c allow for isolated groups
      iso=0
      do n=1,numgroups(i)
       if(itf(n,i).eq.-1)iso=iso+1
      enddo

      ngbig=ng-iso

c     if(ng.gt.nsmall)then
      if(ngbig.gt.nsmall)then
       write(6,*)' There are still more than', nsmall ,'groups' 
       write(6,*)' after maxiters compressions'
       stop
      endif

c     deallocate(ifin)

 
c  compression is finished
      return
      end
