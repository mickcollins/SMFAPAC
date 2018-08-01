      subroutine dbond(x1,x2,r,dr,d2r,nflag)
      implicit double precision(a-h,o-z)

c checked by finite difference 230414

      dimension x1(3),x2(3),dr(2,3),d2r(2,3,2,3),dx(3)

c calculate the bondlength
      r=0.d0
      do k=1,3
       dx(k)=x1(k)-x2(k)
       r=r+dx(k)**2
      enddo
      r2=r
      r=sqrt(r)
      

c first derivatives
      do k=1,3
       dr(1,k)=dx(k)/r
       dr(2,k)=-dr(1,k)
      enddo

c  second derivatives
      if(nflag.eq.2)then

      d2r=0.d0
      do k1=1,3
       d2r(1,k1,1,k1)=1.d0/r-dx(k1)*dr(1,k1)/r2
       d2r(2,k1,2,k1)=1.d0/r+dx(k1)*dr(2,k1)/r2
       d2r(1,k1,2,k1)=-1.d0/r-dx(k1)*dr(2,k1)/r2
       d2r(2,k1,1,k1)=d2r(1,k1,2,k1)
      do k2=1,3
       if(k2.ne.k1)then
        d2r(1,k1,1,k2)=-dx(k1)*dr(1,k2)/r2
        d2r(2,k1,2,k2)=dx(k1)*dr(2,k2)/r2
        d2r(1,k1,2,k2)=-dx(k1)*dr(2,k2)/r2
        d2r(2,k1,1,k2)=dx(k1)*dr(1,k2)/r2
       endif
      enddo
      enddo

      endif

      return
      end
