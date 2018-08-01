      subroutine tsvector(nz,wsave,evec)
      implicit double precision(a-h,o-z)

      integer, allocatable :: nat(:)
      real*8, allocatable :: v(:,:),w(:),work(:),wmax(:),am(:,:)
      dimension wsave(nz),evec(nz+6,nz)

c find which eigenvector has the maximum overlap with the
c atomic motion required

      open(unit=1,file='TSATOMS',status='old')
      read(1,*)
      read(1,*)ntsatoms
      allocate(nat(ntsatoms))
      read(1,*)
      do i=1,ntsatoms
       read(1,*)nat(i)
      enddo
      close(unit=1)
      nts3=3*ntsatoms 
      lwork=max(1,3*nts3-1)
      allocate(v(nts3,nts3))
      allocate(w(nts3))
      allocate(work(lwork))
      allocate(wmax(nz))
c if there are negative eigenvalues, only look at the assoc eigevectors
c else look at them all

      nneg=0
      do i=1,nz
       if(wsave(i).lt.0.d0)nneg=nneg+1
      enddo

      if(nneg.gt.1)then
        go to 100
      else
       go to 200
      endif

100   continue
      write(67,*)' Finding TS vectors from neg eigs'
      do i=1,nz
      wmax(i)=0.d0
      enddo
      do i=1,nz
       if(wsave(i).ge.0.d0)go to 101
       do n1=1,ntsatoms
       do k1=1,3
        j1=3*(n1-1)+k1
        m1=3*(nat(n1)-1)+k1
       do n2=1,ntsatoms
       do k2=1,3
        j2=3*(n2-1)+k2
        m2=3*(nat(n2)-1)+k2
        v(j1,j2)=evec(m1,i)*evec(m2,i)
       enddo
       enddo
       enddo
       enddo
       call dsyev('N','U',nts3,v,nts3,w,work,lwork,ifail)
       al=-1.d6
       do j=1,nts3
        if(w(j).gt.al)al=w(j)
       enddo
       wmax(i)=al
101   enddo
      kts=1
      wts=wmax(1)
      do i=2,nz
       if(wmax(i).gt.wts)then
        wts=wmax(i)
        kts=i
       endif
      enddo
      write(67,*)' The TS vector is number ',kts
      if(kts.ne.1)then
       wstore=wsave(1)
       wsave(1)=wsave(kts)
       wsave(kts)=wstore
       do j=1,nz+6
        estore=evec(j,1)
        evec(j,1)=evec(j,kts)
        evec(j,kts)=estore
       enddo
      endif
      go to 300

200   continue
      write(67,*)' Finding TS vector from all pos eigs'
      do i=1,nz
      wmax(i)=0.d0
      enddo
      do i=1,nz
       do n1=1,ntsatoms
       do k1=1,3
        j1=3*(n1-1)+k1
        m1=3*(nat(n1)-1)+k1
       do n2=1,ntsatoms
       do k2=1,3
        j2=3*(n2-1)+k2
        m2=3*(nat(n2)-1)+k2
        v(j1,j2)=evec(m1,i)*evec(m2,i)
       enddo
       enddo
       enddo
       enddo
       call dsyev('N','U',nts3,v,nts3,w,work,lwork,ifail)
       al=-1.d6
       do j=1,nts3
        if(w(j).gt.al)al=w(j)
       enddo
       wmax(i)=al
      enddo
      kts=1
      wts=wmax(1)
      do i=2,nz
       if(wmax(i).gt.wts)then
        wts=wmax(i)
        kts=i
       endif
      enddo
      write(67,*)' The TS vector is number ',kts
      if(kts.ne.1)then
       wstore=wsave(1)
       wsave(1)=wsave(kts)
       wsave(kts)=wstore
       do j=1,nz+6
        estore=evec(j,1)
        evec(j,1)=evec(j,kts)
        evec(j,kts)=estore
       enddo
      endif

300   continue

c the eigenvector with the greatest overlap with the
c requested atoms is now the first eigenvector

      return
      end
    
