      subroutine hesaugTS(n,fcm,grad1,dgeom1,radius,energy,istep)

      implicit double precision (a-h,o-z)
      dimension grad1(3,n/3),dgeom1(3,n/3)
      dimension fcm(n,n),grad(n),dgeom(n),sol(n-6)
      dimension a(n+1,n+1),v(n+1,n+1),w(n+1),w2(n+1)
      dimension work(10*n)

      dimension gsave(n-6),wsave(n-6),evec(n,n-6)

      character*20 ca,ca1

      parameter (small=2d-5)

      do k=1,3
      do i=1,n/3
       j=3*(i-1)+k
       grad(j)=grad1(k,i)
      enddo
      enddo

      lwork=10*n

c     there will be 6 eigenvalues very close to zero since the 
c       translations and rotations are projected out
c      however there may be some noise and do not want to
c      search along very flat direction, so test for eigenvalues
c      but DO want to search along hydrogen bonds so this
c      threshold is set to 0.0001 which is approx 25 wavenumbers    

c if appropriate get the last predicted energy change
      if(istep.gt.0)then
       open(unit=8,file='predictedenergychange',status='old')
       read(8,*)predicted
       close(unit=8)
       laststep=istep-1
       call filelabel(laststep,ca1)
       n1=index(ca1,' ')-1
       ca='energy.'//ca1(1:n1)
       open(unit=8,file=ca,status='old')
       read(8,*)oldenergy
       close(unit=8)

       actualdeltaE=energy-oldenergy
       write(6,*)' Previously predicted energy change = ',predicted
       write(6,*)' Actual energy change               = ',actualdeltaE
       write(6,*)' The trust radius now = ',radius
      endif


      do j=1,n
      do i=1,n
      v(i,j)=fcm(i,j)
      enddo
      enddo

      call dsyev('V','U',n,v,n+1,w,work,lwork,ifail)

      write(67,*) 'eigenvalues of projected hessian'
      write(67,10) (w(i),i=1,n)
10    format(5f15.5)


c     count the number of non-zero values
      nz=0
      nzero=0
      do i=1,n
      if(dabs(w(i)).le.small)then
      w(i)=0d0
      nzero=nzero+1
      endif
      enddo
      write(67,*) 'There are ', nzero  , ' zero eigenvalues'
      nz=n-nzero

      nneg=0
      do i=1,n
      if(w(i).lt.-small) nneg=nneg+1
      enddo
      write(67,*) 'The hessian has ', nneg ,' negative eigenvalues'
      write(6,*) 'The hessian has ', nneg ,' negative eigenvalues'

c find the simple solution
      m=1
      do i=1,n
       if(abs(w(i)).gt.small)then
        wsave(m)=w(i)
        do j=1,n
         evec(j,m)=v(j,i)
        enddo
        m=m+1
       endif
      enddo

      gnorm=0.d0
      do m=1,nz
       gsave(m)=0.d0
       do i=1,n
        gsave(m)=gsave(m)+evec(i,m)*grad(i)
       enddo
       gnorm=gnorm+gsave(m)**2
      enddo
      gnorm=sqrt(gnorm)

      if(nneg.eq.1)then

       xnorm=0.d0
       do i=1,n
        dgeom(i)=0d0
        do m=1,nz
         dgeom(i)=dgeom(i)-gsave(m)/wsave(m)
        enddo
        xnorm=xnorm+dgeom(i)**2
       enddo
       xnorm=sqrt(xnorm)
      
       if(xnorm.le.radius)then
       write(67,*)' xnorm lt radius and nneg=1'
c  the simple solution is OK, we are done
        go to 500
       endif
      endif

c else we have to follow Bofill 1994

      if(nneg.ne.1)call tsvector(nz,wsave,evec)

c find the shift parameter
      wmin=1.d6
      do i=2,nz
       if(wsave(i).lt.wmin)wmin=wsave(i)
      enddo

600   continue
c initial guess
      uextra=max(wsave(1),-wmin)
      u=gnorm/radius + uextra
      write(67,*)' initial u guess = ',u
      write(67,*)' initial radius = ',radius
c iterate
      do k=1,5
       del=gsave(1)**2/(wsave(1)-u)**2
       ddel=gsave(1)**2/(wsave(1)-u)**3
       do i=2,nz
        del=del+gsave(i)**2/(wsave(i)+u)**2
        ddel=ddel-gsave(i)**2/(wsave(i)+u)**3
       enddo
       del=sqrt(del)
       ddel=ddel/del
       u=u+(1.d0-del/radius)*del/ddel
       write(67,*)' u iterated ', k, u
      enddo
      write(67,*)' Final u = ', u
      if(u.le.uextra)then
       radius=radius/4.d0
       write(67,*)' radius reduced, calc u again'
       write(67,*)' radius now = ',radius
       go to 600
      endif
c evaluate the solutions for the coefficients
      write(67,*)' Solution for coefficients'
      sol(1)=-gsave(1)/(wsave(1)-u)
      write(67,*)' s ',1,1.d0/(wsave(1)-u),gsave(1)
      do i=2,nz
       sol(i)=-gsave(i)/(wsave(i)+u)
       write(67,*)' s ',i,1.d0/(wsave(i)+u),gsave(i)
      enddo
c evaluate dgeom
      do i=1,n
       dgeom(i)=sol(1)*evec(i,1)
      enddo
      do i=1,n
       do m=2,nz
        dgeom(i)=dgeom(i)+sol(m)*evec(i,m)
       enddo
      enddo
       xnorm=0.d0
       do i=1,n
        xnorm=xnorm+dgeom(i)**2
       enddo
       xnorm=sqrt(xnorm)
       write(67,*)' xnorm, radius ',xnorm, radius

500   continue

c     scale using trust radius
c     NB the trust radius should really be set dynamically
c simpler than Anglada and Bofill
      if(istep.gt.1)then
       ratio=actualdeltaE/predicted
       if(xnorm.gt.0.95d0*radius)then
        if(ratio.gt.0.85.and.ratio.lt.1.15)then
         radius=radius*1.4142d0
        else
         if(ratio.gt.0.75.and.ratio.lt.1.25)radius=radius/2.d0
        endif
       endif

       open(unit=8,file='TRUSTRADIUS',status='old')
        write(8,*)radius
       close(unit=8)
      endif

      write(67,*)' end of hesaugTS'
      write(67,*)

      do k=1,3
      do i=1,n/3
       j=3*(i-1)+k
       dgeom1(k,i)=dgeom(j)
      enddo
      enddo

      return
      end

