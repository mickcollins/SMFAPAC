      program PAA 

c  integrates imag frew polarizability products to get 
c  PAA factor of some molecule

      implicit double precision(a-h,o-z)

      dimension v(16),t(16),wt(16)

      dimension a(3,3,16),as(3,3),aa(3,3,16),asper(16)

      dimension c(3,3)

      character*2 lab(3)

      v0=0.2d0

      pi2=4.d0*acos(0.d0)

      x1=-1.d0
      x2=1.d0

      call gauleg(x1,x2,t,wt,16)

c input file (unit 5) should have
c  number of atoms
c  coordinates for each atom as label x,y,z
c  xx header
c  xx imag freq pols at each freq
c  yx header
c  yx imag freq pols at each freq
c  yy header
c  yy imag freq pols at each freq
c  zx header
c  zx imag freq pols at each freq
c  zy header
c  zy imag freq pols at each freq
c  zz header
c  zz imag freq pols at each freq

c     open(unit=1,file='dispersion_all',status='old')

c     do m=1,13
      do m=1,1

      read(5,*)natom
      write(6,*)natom
      write(6,*)' Coordinates'

      do n=1,natom
       read(5,*)lab(n),(c(n,j),j=1,3)
       write(6,100)lab(n),(c(n,j),j=1,3)
      enddo
100   format(a2,3f13.6)
c100   format(a2,3f19.14)
      do i1=1,3
      do i2=1,i1

       read(5,*)
       read(5,*)aaa,as(i1,i2)
      do k=1,16
       read(5,*)rub,a(i1,i2,k)
      enddo

      enddo
      enddo

c average static polarizability
      asaver=(as(1,1)+as(2,2)+as(3,3))/3.d0      

      do k=1,16
       asper(17-k)=(a(1,1,k)+a(2,2,k)+a(3,3,k))/3.d0
      enddo

c scale out the static value
      
      do k=1,16
       asper(k)=asper(k)/asaver
      enddo

      sum=0.d0
      do k=1,16
       sum=sum+wt(k)*asper(k)*asper(k)/(1.d0+t(k))**2
      enddo

      sum=sum*2.d0*v0/pi2

      write(6,*)' PAA factor'
      write(6,*)sum
c end loop over geometries
      enddo

      end

      SUBROUTINE gauleg(x1,x2,x,w,n)
      implicit double precision(a-h,o-z)
      INTEGER n
c     REAL*8 x1,x2,x(n),w(n)

      dimension x(n),w(n)

      PARAMETER (EPS=3.d-14)
c     INTEGER i,j,m
c     DOUBLE PRECISION p1,p2,p3,pp,xl,xm,z,z1
      m=(n+1)/2
      xm=0.5d0*(x2+x1)
      xl=0.5d0*(x2-x1)
      do 12 i=1,m
        z=cos(3.141592654d0*(i-.25d0)/(n+.5d0))
1       continue
          p1=1.d0
          p2=0.d0
          do 11 j=1,n
            p3=p2
            p2=p1
            p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
11        continue
          pp=n*(z*p1-p2)/(z*z-1.d0)
          z1=z
          z=z1-p1/pp
        if(abs(z-z1).gt.EPS)goto 1
        x(i)=xm-xl*z
        x(n+1-i)=xm+xl*z
        w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
        w(n+1-i)=w(i)
12    continue
      return
      END
