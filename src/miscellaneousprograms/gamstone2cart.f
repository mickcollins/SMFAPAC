      program gamstone2cart
      implicit double precision(a-h,o-z)

      real*8, allocatable :: c(:,:),ch(:),d(:,:),q(:,:,:),oct(:,:,:,:)
      real*8, allocatable :: chZ(:)


      character*2, allocatable :: lab(:)

      dimension q0(6),o0(10),r(3,3,3)

      character*3 rub

      bohr=1.d0/1.8897259886d0

      read(5,*)n
      natom=n/2

      allocate(c(natom,3))
      allocate(ch(natom))
      allocate(d(natom,3))
      allocate(q(natom,3,3))
      allocate(oct(natom,3,3,3))
      allocate(lab(natom))
      allocate(chZ(natom))
      
      do n=1,natom
       read(5,*)i1,rub,chZ(n),(c(n,k),k=1,3)
      enddo
      do n=1,natom
       read(5,*)i1,rub,ch(n)
       ch(n)=ch(n)+chZ(n)
      enddo

c convert coords to Angs
      do n=1,natom
      do k=1,3
       c(n,k)=c(n,k)*bohr
      enddo
      enddo

      do n=1,natom
       read(5,100)lab(n),(d(n,k),k=1,3)
      enddo
100   format(1x,a2,6x,3f14.5)
      do n=1,natom
       read(5,101)lab(n),(q0(k),k=1,6)
c convert
       trce=q0(1)+q0(2)+q0(3)
       q(n,1,1)=0.5d0*(3.d0*q0(1)-trce)
       q(n,2,2)=0.5d0*(3.d0*q0(2)-trce)
       q(n,3,3)=0.5d0*(3.d0*q0(3)-trce)
       q(n,1,2)=1.5d0*q0(4)
        q(n,2,1)=q(n,1,2)
       q(n,1,3)=1.5d0*q0(5)
        q(n,3,1)=q(n,1,3)
       q(n,2,3)=1.5d0*q0(6)
        q(n,3,2)=q(n,2,3)
      enddo
101   format(1x,a2,6x,6f11.5)
      do n=1,natom
       read(5,101)lab(n),(o0(k),k=1,6)
       read(5,102)(o0(k),k=7,10)
       r(1,1,1)=o0(1)
       r(2,2,2)=o0(2)
       r(3,3,3)=o0(3)
       r(1,1,2)=o0(4)
        r(1,2,1)=r(1,1,2)
        r(2,1,1)=r(1,1,2)
       r(1,1,3)=o0(5)
        r(1,3,1)=r(1,1,3)
        r(3,1,1)=r(1,1,3)
       r(1,2,2)=o0(6)
        r(2,1,2)=r(1,2,2)
        r(2,2,1)=r(1,2,2)
       r(2,2,3)=o0(7)
        r(2,3,2)=r(2,2,3)
        r(3,2,2)=r(2,2,3)
       r(1,3,3)=o0(8)
        r(3,1,3)=r(1,3,3)
        r(3,3,1)=r(1,3,3)
       r(2,3,3)=o0(9)
        r(3,2,3)=r(2,3,3)
        r(3,3,2)=r(2,3,3)
       r(1,2,3)=o0(10)
        r(2,1,3)=r(1,2,3)
        r(3,2,1)=r(1,2,3)
        r(1,3,2)=r(1,2,3)
        r(2,3,1)=r(1,2,3)
        r(3,1,2)=r(1,2,3)
      do k1=1,3
      do k2=1,3
      do k3=1,3
       oct(n,k1,k2,k3)=5.d0*r(k1,k2,k3)/2.d0
       if(k2.eq.k3)then
        oct(n,k1,k2,k3)=oct(n,k1,k2,k3)
     .      -0.5d0*(r(1,1,k1)+r(2,2,k1)+r(3,3,k1))
       endif
       if(k1.eq.k3)then
        oct(n,k1,k2,k3)=oct(n,k1,k2,k3)
     .      -0.5d0*(r(1,1,k2)+r(2,2,k2)+r(3,3,k2))
       endif
       if(k1.eq.k2)then
        oct(n,k1,k2,k3)=oct(n,k1,k2,k3)
     .      -0.5d0*(r(1,1,k3)+r(2,2,k3)+r(3,3,k3))
       endif
      enddo
      enddo
      enddo

      enddo
102   format(9x,4f11.5)


c now make a ".cart" file
      write(6,*)"  The distributed Cartesian multipoles"
      write(6,*)natom
      do n=1,natom
       write(6,201)lab(n)
201    format(a2)
       do k=1,3
        write(6,*)c(n,k)
       enddo
      write(6,*)"  Charge"
      write(6,*)ch(n)
      write(6,*)"  Dipole"
       do k=1,3
        write(6,*)d(n,k)
       enddo
      write(6,*)"  Quodrupole"
      do k1=1,3
      do k2=k1,3
       write(6,*)q(n,k1,k2)
      enddo
      enddo
      write(6,*)"  Octapole"
      do k1=1,3
      do k2=k1,3
      do k3=k2,3
       write(6,*)oct(n,k1,k2,k3)
      enddo
      enddo
      enddo
      write(6,*)"  Hexadecapole"
      do k1=1,3
      do k2=k1,3
      do k3=k2,3
      do k4=k3,3
       write(6,*)0.d0
      enddo
      enddo
      enddo
      enddo

      enddo


      end
