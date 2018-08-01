      subroutine fragsteps

      use fractdata
      implicit double precision(a-h,o-z)

      integer, allocatable  :: ilinkstore(:,:)
      allocate(ilinkstore(natom,natom))
      allocate(ilink(nsmall,nsmall))
      allocate(ma(natom,natom))
 
10    continue
      it=0
      do n=1,numat(nfrag)
      if(itype(nfrag,n).gt.it)it=itype(nfrag,n)
      enddo
      if(numat(nfrag).le.Level+it)nstop(nfrag)=1

      if(nstop(nfrag).gt.0)go to 1
      call shrink
      do n1=1,numat(nfrag)
      do n2=1,numat(nfrag)
       ilinkstore(n1,n2)=ilink(n1,n2)
      enddo
      enddo


      call breakgroup


      deallocate(mats)
      deallocate(matr)
      deallocate(itp)
      deallocate(ibo)


      if(nstop(nfrag).gt.0)go to 1
      allocate(notthere(natom))
c first eliminate nelim
      do n=1,numat(nfrag)
       ilink(nelim,n)=0
       ilink(n,nelim)=0
       notthere(n)=1
      enddo
      notthere(nelim)=0
      call findfragelim
      isgn=1
      call makefrag(isgn)

c then eliminate the far groups

      do n1=1,numat(nfrag)
      do n2=1,numat(nfrag)
       ilink(n1,n2)=ilinkstore(n1,n2)
      enddo
      enddo

      do n1=1,numat(nfrag)
       if(ma(nelim,n1).eq.0)then
        do n2=1,numat(nfrag)
         ilink(n1,n2)=0
         ilink(n2,n1)=0
        enddo
        notthere(n1)=0
       endif
      enddo
      notthere(nelim)=1
      call findfragelim
      isgn=1
      call makefrag(isgn)

c now eliminate both
      do n=1,numat(nfrag)
       ilink(nelim,n)=0
       ilink(n,nelim)=0
      enddo
      notthere(nelim)=0
      call findfragelim
      isgn=-1
      call makefrag(isgn)
      nstop(nfrag)=2
      deallocate(notthere)

1      continue

c     if(nf.gt.(nfragm-3*nsmall))then
      if(nf.gt.(nfragm-natom))then
        write(6,*)'too many fragments'
      write(6,*)nfrag,nf
      call cancel
      write(6,*)nfrag,nf
c     if(nfrag.eq.1.and.nf.gt.(nfragm-3*nsmall))then
      if(nfrag.eq.1.and.nf.gt.(nfragm-natom))then
       write(6,*)' Fragmentation will freeze, so abort'
       stop
      endif
      else
      nfrag=nfrag+1
      endif

      if(nfrag.gt.nf)go to 2

c     stop
      go to 10
2     continue

      call cancel

      deallocate(ilinkstore)
      deallocate(ilink)
      deallocate(ma)

      return
      end
