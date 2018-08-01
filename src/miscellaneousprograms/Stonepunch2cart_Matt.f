      program Stonepunch2c

c  this program reads in a "molecule.punch" file
c  as produced by gdm2.2.2.03 to get the distributed
c  multipole moments
c  The spherical tensor data is converted to
c  Cartesian tensor data

      implicit double precision(a-h,o-z)
      parameter(natomm=500,maxrk=5)
      parameter(CVS=1.8897259886d0)

      character*2 lab(natomm)
      character*6 com1
      character*8 com2

      dimension c(natomm,3)
      dimension irk(natomm)

c  we assume maximum rank4 tensors

      dimension q(natomm,0:maxrk,2*maxrk+1)

      dimension ch(natomm),d(natomm,3),qu(natomm,3,3),
     .          oct(natomm,3,3,3),hex(natomm,3,3,3,3)

      do n=1,natomm
      ch(n)=0.d0
      do k1=1,3
      d(n,k1)=0.d0
      do k2=1,3
      qu(n,k1,k2)=0.d0
      do k3=1,3
      oct(n,k1,k2,k3)=0.d0
      do k4=1,3
      hex(n,k1,k2,k3,k4)=0.d0
      enddo
      enddo
      enddo
      enddo
      enddo

      read(5,*)
      read(5,*)

      n=1
1     continue

      read(5,*,end=2)
      read(5,100,end=2)lab(n),com1,(c(n,k),k=1,3),com2,irank
      irk(n)=irank
      do m=0,irank
      m2=2*m+1
c     MAA Edit 14/03/09 changed read(5,*) to read(5,'(5(f16.10))') 
      read(5,*)(Q(n,m,k),k=1,m2)
c     read(5,'(5(f16.10))')(Q(n,m,k),k=1,m2)
c     write(6,*)(Q(n,m,k),k=1,m2)
      enddo
c     stop
c  convert to Cartesian data

      m=0
      ch(n)=Q(n,0,1)
      if(irank.eq.m)go to 3
      m=1
      d(n,1)=Q(n,1,2)
      d(n,2)=Q(n,1,3)
      d(n,3)=Q(n,1,1)
      if(irank.eq.m)go to 3
      m=2
      qu(n,1,1)=0.5d0*(sqrt(3.d0)*Q(n,2,4)-Q(n,2,1))
      qu(n,2,2)=-0.5d0*(sqrt(3.d0)*Q(n,2,4)+Q(n,2,1))
      qu(n,3,3)=Q(n,2,1)
      qu(n,1,2)=0.5d0*sqrt(3.d0)*q(n,2,5)
      qu(n,2,1)=qu(n,1,2)
      qu(n,1,3)=0.5d0*sqrt(3.d0)*q(n,2,2)
      qu(n,3,1)=qu(n,1,3)
      qu(n,2,3)=0.5d0*sqrt(3.d0)*q(n,2,3)
      qu(n,3,2)=qu(n,2,3)
      if(irank.eq.m)go to 3
      m=3
      s58=sqrt(5.d0/8.d0)
      s38=sqrt(3.d0/8.d0)
      s12=1.d0/sqrt(24.d0)
      s512=sqrt(5.d0/12.d0)
      s23=sqrt(2.d0/3.d0)

      oct(n,1,1,1)=s58*Q(n,3,6)-s38*Q(n,3,2)
      oct(n,1,1,2)=s58*Q(n,3,7)-s12*Q(n,3,3)
      oct(n,1,2,2)=-s58*Q(n,3,6)-s12*Q(n,3,2)
      oct(n,2,2,2)=-s58*Q(n,3,7)-s38*Q(n,3,3)
      oct(n,1,1,3)=s512*Q(n,3,4)-0.5d0*Q(n,3,1)
      oct(n,1,2,3)=s512*Q(n,3,5)
      oct(n,2,2,3)=-s512*Q(n,3,4)-0.5d0*Q(n,3,1)
      oct(n,1,3,3)=s23*Q(n,3,2)
      oct(n,2,3,3)=s23*Q(n,3,3)
      oct(n,3,3,3)=Q(n,3,1)

c  symmetrize
c     do k1=1,3
c     do k2=1,k1
c     do k3=1,k2
c      oct(n,k1,k2,k3)=oct(n,k3,k2,k1)
c      oct(n,k1,k3,k2)=oct(n,k3,k2,k1)
c      oct(n,k2,k1,k3)=oct(n,k3,k2,k1)
c      oct(n,k2,k3,k1)=oct(n,k3,k2,k1)
c      oct(n,k3,k1,k2)=oct(n,k3,k2,k1)
c     enddo
c     enddo
c     enddo

      if(irank.eq.m)go to 3
      m=4

      f38=3.d0/8.d0
      f4=0.25d0
      f32=1.d0/32.d0
      f8=1.d0/8.d0
      f16=1.d0/16.d0
      s5=sqrt(5.d0)
      s35=sqrt(35.d0)
      s10=sqrt(10.d0)
      s70=sqrt(70.d0)
      s58=sqrt(5.d0/8.d0)

      hex(n,1,1,1,1)=f38*Q(n,4,1)-f4*s5*Q(n,4,4)+f8*s35*Q(n,4,8)
      hex(n,1,1,1,2)=f8*(-s5*Q(n,4,5)+s35*Q(n,4,9))
      hex(n,1,1,2,2)=f8*(Q(n,4,1)-s35*Q(n,4,8))
      hex(n,1,2,2,2)=-f8*(s5*Q(n,4,5)+s35*Q(n,4,9))
      hex(n,2,2,2,2)=f38*Q(n,4,1)+f4*s5*Q(n,4,4)+f8*s35*Q(n,4,8)
      hex(n,1,1,1,3)=f16*(-3.d0*s10*Q(n,4,2)+s70*Q(n,4,6))
      hex(n,1,1,2,3)=f16*(-s10*Q(n,4,3)+s70*Q(n,4,7))
      hex(n,1,2,2,3)=-f16*(s10*Q(n,4,2)+s70*Q(n,4,6))
      hex(n,2,2,2,3)=-f16*(3*s10*Q(n,4,3)+s70*Q(n,4,7))
      hex(n,1,1,3,3)=-0.5d0*Q(n,4,1)+f4*s5*Q(n,4,4)
      hex(n,1,2,3,3)=f4*s5*Q(n,4,5)
      hex(n,2,2,3,3)=-0.5d0*Q(n,4,1)-f4*s5*Q(n,4,4)
      hex(n,1,3,3,3)=s58*Q(n,4,2)
      hex(n,2,3,3,3)=s58*Q(n,4,3)
      hex(n,3,3,3,3)=Q(n,4,1)

c  symmetrize
c     do k1=1,3
c     do k2=1,k1
c     do k3=1,k2
c     do k4=1,k3
c     hex(n,k4,k3,k1,k2)=hex(n,k4,k3,k2,k1)
c     hex(n,k4,k2,k1,k3)=hex(n,k4,k3,k2,k1)
c     hex(n,k4,k2,k3,k1)=hex(n,k4,k3,k2,k1)
c     hex(n,k4,k1,k2,k3)=hex(n,k4,k3,k2,k1)
c     hex(n,k4,k1,k3,k2)=hex(n,k4,k3,k2,k1)

c     hex(n,k3,k1,k2,k4)=hex(n,k4,k3,k2,k1)
c     hex(n,k3,k1,k4,k2)=hex(n,k4,k3,k2,k1)
c     hex(n,k3,k2,k1,k4)=hex(n,k4,k3,k2,k1)
c     hex(n,k3,k2,k4,k1)=hex(n,k4,k3,k2,k1)
c     hex(n,k3,k4,k2,k1)=hex(n,k4,k3,k2,k1)
c     hex(n,k3,k4,k1,k2)=hex(n,k4,k3,k2,k1)

c     hex(n,k2,k1,k3,k4)=hex(n,k4,k3,k2,k1)
c     hex(n,k2,k1,k4,k3)=hex(n,k4,k3,k2,k1)
c     hex(n,k2,k3,k1,k4)=hex(n,k4,k3,k2,k1)
c     hex(n,k2,k3,k4,k1)=hex(n,k4,k3,k2,k1)
c     hex(n,k2,k4,k1,k3)=hex(n,k4,k3,k2,k1)
c     hex(n,k2,k4,k3,k1)=hex(n,k4,k3,k2,k1)

c     hex(n,k1,k2,k3,k4)=hex(n,k4,k3,k2,k1)
c     hex(n,k1,k2,k4,k3)=hex(n,k4,k3,k2,k1)
c     hex(n,k1,k3,k2,k4)=hex(n,k4,k3,k2,k1)
c     hex(n,k1,k3,k4,k2)=hex(n,k4,k3,k2,k1)
c     hex(n,k1,k4,k2,k3)=hex(n,k4,k3,k2,k1)
c     hex(n,k1,k4,k3,k2)=hex(n,k4,k3,k2,k1)

c     enddo
c     enddo
c     enddo
c     enddo

3     continue
      n=n+1
      go to 1
2     continue

      natom=n-1

      write(6,*)' The distributed Cartesian multipoles'
      write(6,*)natom
      do n=1,natom
c     write(6,*)' Coordinates ',n
c     write(6,*)(c(n,k)*CVS,k=1,3)
      write(6,200)lab(n)
200   format(a2)
c     write(6,*)irk(n)
      do k=1,3
c     write(6,*)c(n,k)*CVS
      write(6,*)c(n,k)
      enddo
      write(6,*)' Charge'
      write(6,*)ch(n)
      !enddo ! MAA Added
c     if(irk(n).eq.0)go to 10
      write(6,*)' Dipole'
      do k1=1,3
      write(6,*)d(n,k1)
      enddo
c     if(irk(n).eq.1)go to 10
      write(6,*)' Quodrupole'
      do k1=1,3
      do k2=k1,3
      write(6,*)qu(n,k1,k2)
      enddo
      enddo
c     if(irk(n).eq.2)go to 10


      write(6,*)' Octapole'
      do k1=1,3
      do k2=k1,3
      do k3=k2,3
      write(6,*)oct(n,k1,k2,k3)
      enddo
      enddo
      enddo
c     if(irk(n).eq.3)go to 10
      write(6,*)' Hexadecapole'
      do k1=1,3
      do k2=k1,3
      do k3=k2,3
      do k4=k3,3
      write(6,*)hex(n,k1,k2,k3,k4)
      enddo
      enddo
      enddo
      enddo

10    continue

      enddo

c     write(6,*)' Hexadecapole'
c     do k1=1,3
c     do k2=k1,3
c     do k3=k2,3
c     do k4=k3,3
c     write(6,400)k1,k2,k3,k4,hex(1,k1,k2,k3,k4)
c     enddo
c     enddo
c     enddo
c     enddo
400   format(1x,4i3,e15.6)

c     write(6,*)' The distributed Cartesian multipoles'
c     write(6,*)natom
c     do n=1,natom
c     write(6,101)lab(n),(c(n,k)*CVS,k=1,3)
c     write(6,*)irk(n)
c     if(irk(n).ge.0)write(6,*)ch(n)
c     if(irk(n).ge.1)write(6,*)(d(n,k),k=1,3)
c     if(irk(n).ge.2)write(6,*)((qu(n,k1,k2),k2=k1,3),k1=1,3)
c     if(irk(n).ge.3)then
c     write(6,*)(((oct(n,k1,k2,k3),k3=k2,3),k2=k1,3),k1=1,3)
c     endif
c     if(irk(n).ge.4)then
c     write(6,*)
c    . ((((hex(n,k1,k2,k3,k4),k4=k3,3),k3=k2,3),k2=k1,3),k1=1,3)
c     endif
c     enddo


100   format(a2,a6,3F16.10,a8,i3)
101   format(a2,3F16.10)

      end


