      program newHatoms
      implicit double precision(a-h,o-z)

      real*8, allocatable :: c(:,:)
      character*2, allocatable :: lab(:)
      character*80 filename

      dimension x(3),y(3),z(3),h(3)

      open(unit=1,file='xyzFILENAME',status='old')
      read(1,100)filename
      close(unit=1)
100   format(a80)
      open(unit=1,file=trim(filename),status='old')
      read(1,*)natom
      allocate(lab(natom))
      allocate(c(natom,3))
      read(1,*)
      do n=1,natom
       read(1,*)lab(n),(c(n,k),k=1,3)
      enddo
      close(unit=1)

      read(5,*) n1,n2,n3

      sumx=0.d0
      do k=1,3
       x(k)=c(n1,k)-c(n2,k)
       y(k)=c(n3,k)-c(n2,k)
       sumx=sumx+x(k)**2
      enddo
      sumx=sqrt(sumx)
      do k=1,3
      x(k)=x(k)/sumx
      enddo
      dot=0.d0
      do k=1,3
       dot=dot+x(k)*y(k)
      enddo
      sumy=0.d0
      do k=1,3
       y(k)=y(k)-dot*x(k)
       sumy=sumy+y(k)**2
      enddo
      sumy=sqrt(sumy)
      do k=1,3
       y(k)=y(k)/sumy
      enddo


      z(1)=x(2)*y(3)-x(3)*y(2)
      z(2)=x(3)*y(1)-x(1)*y(3)
      z(3)=x(1)*y(2)-x(2)*y(1)

c     sumx=0.d0
c     sumy=0.d0
c     dot=0.d0
c     do k=1,3
c      sumx=sumx+x(k)**2
c      sumy=sumy+y(k)**2
c      dot=dot+x(k)*y(k)
c     enddo
c     dz=0.d0
c     dot1=0.d0
c     dot2=0.d0
c     do k=1,3
c      dz=dz+z(k)**2
c      dot1=dot1+x(k)*z(k)
c      dot2=dot2+y(k)*z(k)
c     enddo
c     write(6,*)sumx,sumy,dz
c     write(6,*)dot,dot1,dot2

      read(5,*) theta, psi

      pi=2.d0*acos(0.d0)
      theta=theta*pi/180.0
      psi=psi*pi/180.d0
      ct=cos(theta)
      st=sqrt(1.d0-ct**2)
      cp=cos(psi)
      sp=sin(psi)
      do k=1,3
       h(k)=ct*z(k)+st*cp*x(k)+st*sp*y(k)+c(n2,k)
      enddo

c     sumh=0.d0
c      do k=1,3
c       sumh=sumh+h(k)**2
c      enddo
c      write(6,*)sumh

      write(6,200)'H ',(h(k),k=1,3)
200   format(a2,3f13.6)

      end
