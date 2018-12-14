      subroutine extrafrags
      use fractdata
      implicit double precision(a-h,o-z)
      write(6,*)'calling makeneighbours'
      write(6,*)' nf = ',nf
      write(6,*)' nfragm = ',nfragm
      write(6,*)' nbondsextra = ',nbondsextra
      call makeneighbours
      write(6,*)' finished makeneighbours'

      nfrag=1
c      allocate(map(nsmall,nsmall))
c      allocate(kf(nsmall))
      call fragsteps
      write(6,*)' finished fragsteps'
      call finalcancel
      write(6,*)' finished finalcancel'

      return
      end


      subroutine makeneighbours
      use fractdata
      implicit double precision(a-h,o-z)

      nsmallkat=10*nsmall
      allocate(kat(nsmallkat))

      do k=1,nbondsextra
       m1=mbn(k)
       m2=nbn(k)
       kat(1)=m1
       kat(2)=m2
       nk=2

       do i=1,Level
        ik=nk
        do n=1,nk
         if(itf(kat(n),1).gt.-1)then
          do m=1,itf(kat(n),1)+1
           ik=ik+1
           if(ik.gt.nsmallkat)then
            write(6,*)' ik is too big'
            stop
           endif
           kat(ik)=iabs(ibf(kat(n),m,1))
          enddo
         endif
        enddo
        nk=ik
       enddo
c cancel repeats
       do j1=1,nk-1
        if(kat(j1).eq.0)go to 2
       do j2=j1+1,nk
        if(kat(j2).eq.0)go to 1
        if(kat(j1)-kat(j2).eq.0)kat(j2)=0
1      enddo
2      enddo
       ic=0
       do j1=1,nk
        if(kat(j1).ne.0)then
         ic=ic+1
         kat(ic)=kat(j1)
        endif
       enddo

       nk=ic
       write(81,*)(kat(kk),kk=1,nk)
       write(81,*)
       if(nk.gt.nsmall)then
        write(6,*)' nk.gt.nsmall, nk = ',nk
        write(6,*)' at extrabond number ',k
        write(6,*)' for atoms '
        write(6,*)nfam(m1)+nfam(m2)
        write(6,*)nfam(m1),nfam(m2)
        do i1=1,nfam(m1)
        j1=ifam(m1,i1)
        write(6,100)atoms(j1),(c(j1,l),l=1,3)
        enddo
        do i1=1,nfam(m2)
        j1=ifam(m2,i1)
        write(6,100)atoms(j1),(c(j1,l),l=1,3)
        enddo
        stop
       endif
100    format(a2,3f13.4)

       call makeextrafrags(nk,m1,m2,multn(k))

c close k loop
      enddo

      return
      end

      subroutine makeextrafrags(nk,m1,m2,ksign)
      use fractdata
      implicit double precision(a-h,o-z)

      do n=1,nsmall
      do m=1,nsmall
       ibond(n,m)=0
      enddo
      enddo

c subtract the non-bonded
      nf=nf+1
      if(nf.ge.nfragm)then
       write(6,*)' nf.ge.nfragm in makeextrafrags'
       stop
      endif
      nstop(nf)=0
      isign(nf)=-1
      numat(nf)=nk
      do i=1,nk
       natstore(nf,i)=kat(i)
      enddo
      do n=1,numat(nf)
       itype(nf,n)=0
       if(itf(kat(n),1).gt.-1)then
        do i=1,itf(kat(n),1)+1
        do m=1,numat(nf)
         if(iabs(ibf(kat(n),i,1)).eq.kat(m))then
          itype(nf,n)=itype(nf,n)+1
         ibond(n,itype(nf,n))=m*ibf(kat(n),i,1)/iabs(ibf(kat(n),i,1))
         endif
        enddo
        enddo
       endif
       itype(nf,n)=itype(nf,n)-1
      enddo

      ic=0
      do n=1,numat(nf)
      do m=1,itype(nf,n)+1
       ic=ic+1
       ib1(nf,ic)=ibond(n,m)
      enddo
      enddo

      write(81,*)' nf = ',nf
      do n=1,numat(nf)
      write(81,*)' natstore = ',natstore(nf,n)
      write(81,*)(ibond(n,m),m=1,itype(nf,n)+1)
      enddo
      write(81,*)
      write(81,*)

c update itf and ibf
      itf(m1,1)=itf(m1,1)+1
      ibf(m1,itf(m1,1)+1,1)=ksign*m2
      itf(m2,1)=itf(m2,1)+1
      ibf(m2,itf(m2,1)+1,1)=ksign*m1

      do n=1,nsmall
      do m=1,nsmall
       ibond(n,m)=0
      enddo
      enddo

c repeat with this extra bond
      nf=nf+1
      if(nf.ge.nfragm)then
       write(6,*)' nf.ge.nfragm in makeextrafrags'
       stop
      endif
      nstop(nf)=0
      isign(nf)=1
      numat(nf)=nk
      do i=1,nk
       natstore(nf,i)=kat(i)
      enddo
      do n=1,numat(nf)
       itype(nf,n)=0
       if(itf(kat(n),1).gt.-1)then
        do i=1,itf(kat(n),1)+1
        do m=1,numat(nf)
         if(iabs(ibf(kat(n),i,1)).eq.kat(m))then
          itype(nf,n)=itype(nf,n)+1
         ibond(n,itype(nf,n))=m*ibf(kat(n),i,1)/iabs(ibf(kat(n),i,1))
         endif
        enddo
        enddo
       endif
      itype(nf,n)=itype(nf,n)-1
      enddo

      ic=0
      do n=1,numat(nf)
      do m=1,itype(nf,n)+1
       ic=ic+1
       ib1(nf,ic)=ibond(n,m)
      enddo
      enddo

      write(81,*)' nf = ',nf
      do n=1,numat(nf)
      write(81,*)' natstore = ',natstore(nf,n)
      write(81,*)(ibond(n,m),m=1,itype(nf,n)+1)
      enddo
      write(81,*)
      write(81,*)
      
      return
      end

