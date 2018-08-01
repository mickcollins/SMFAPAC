      subroutine dangle(x1,x2,x3,a,da,d2a,nflag)
      implicit double precision(a-h,o-z)

c checked by finite difference 230414

      dimension x1(3),x2(3),x3(3),da(3,3),d2a(3,3,3,3)
      dimension dr1(2,3),dr2(2,3),dr3(2,3)
      dimension d2r1(2,3,2,3),d2r2(2,3,2,3),d2r3(2,3,2,3)

      dimension df1(3,3),df2(3,3),df3(3,3)
      dimension d2f1(3,3,3,3),d2f2(3,3,3,3),d2f3(3,3,3,3)

c the angle x1...x2...x3

c actually calculate cos of the angle
c and derivatives of cos

      call dbond(x1,x2,r1,dr1,d2r1,nflag)
      call dbond(x3,x2,r2,dr2,d2r2,nflag)
      call dbond(x1,x3,r3,dr3,d2r3,nflag)


      denom=2.d0*r1*r2
      denom2=denom**2
      anom=r1**2+r2**2-r3**2

c cos of the angle
      a=anom/denom

      do k=1,3
       df1(1,k)=dr1(1,k)
       df1(2,k)=dr1(2,k)
       df1(3,k)=0.d0
       df2(1,k)=0.d0
       df2(2,k)=dr2(2,k)
       df2(3,k)=dr2(1,k)
       df3(1,k)=dr3(1,k)
       df3(2,k)=0.d0
       df3(3,k)=dr3(2,k)
      enddo

      do k1=1,3
      do k2=1,3

       d2f1(1,k1,1,k2)=d2r1(1,k1,1,k2)
       d2f1(1,k1,2,k2)=d2r1(1,k1,2,k2)
       d2f1(2,k1,1,k2)=d2r1(1,k1,2,k2)
       d2f1(1,k1,3,k2)=0.d0
       d2f1(3,k1,1,k2)=0.d0
       d2f1(2,k1,2,k2)=d2r1(2,k1,2,k2)
       d2f1(2,k1,3,k2)=0.d0
       d2f1(3,k1,2,k2)=0.d0
       d2f1(3,k1,3,k2)=0.d0

       d2f2(1,k1,1,k2)=0.d0
       d2f2(1,k1,2,k2)=0.d0
       d2f2(2,k1,1,k2)=0.d0
       d2f2(1,k1,3,k2)=0.d0
       d2f2(3,k1,1,k2)=0.d0
       d2f2(2,k1,2,k2)=d2r2(2,k1,2,k2)
       d2f2(2,k1,3,k2)=d2r2(1,k1,2,k2)
       d2f2(3,k1,2,k2)=d2r2(1,k1,2,k2)
       d2f2(3,k1,3,k2)=d2r2(1,k1,1,k2)

       d2f3(1,k1,1,k2)=d2r3(1,k1,1,k2)
       d2f3(1,k1,2,k2)=0.d0
       d2f3(2,k1,1,k2)=0.d0
       d2f3(1,k1,3,k2)=d2r3(1,k1,2,k2)
       d2f3(3,k1,1,k2)=d2r3(1,k1,2,k2)
       d2f3(2,k1,2,k2)=0.d0
       d2f3(2,k1,3,k2)=0.d0
       d2f3(3,k1,2,k2)=0.d0
       d2f3(3,k1,3,k2)=d2r3(2,k1,2,k2)

      enddo
      enddo



c first derivatives
      do k=1,3
      do n=1,3
       da(n,k)=(r1*df1(n,k)+r2*df2(n,k)-r3*df3(n,k)
     .         -df1(n,k)*r2*a-r1*df2(n,k)*a)/(r1*r2)
      enddo
      enddo

c second derivatives
      if(nflag.eq.2)then

      do k1=1,3
      do k2=1,3
      do n1=1,3
      do n2=1,3

      d2a(n1,k1,n2,k2)=(df1(n1,k1)*df1(n2,k2)+r1*d2f1(n1,k1,n2,k2)
     .                 +df2(n1,k1)*df2(n2,k2)+r2*d2f2(n1,k1,n2,k2)
     .                 -df3(n1,k1)*df3(n2,k2)-r3*d2f3(n1,k1,n2,k2)
     .                 -d2f1(n1,k1,n2,k2)*r2*a-df1(n1,k1)*df2(n2,k2)*a
     .                 -df1(n1,k1)*r2*da(n2,k2)
     .                 -d2f2(n1,k1,n2,k2)*r1*a-df2(n1,k1)*df1(n2,k2)*a
     .                 -df2(n1,k1)*r1*da(n2,k2)
     .                 -df1(n2,k2)*r2*da(n1,k1)-df2(n2,k2)*r1*da(n1,k1))
     .                 /(r1*r2)

      enddo
      enddo
      enddo
      enddo

      endif

      return
      end
