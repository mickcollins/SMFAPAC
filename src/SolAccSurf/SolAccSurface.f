      program SolAccSurface

      implicit double precision(a-h,o-z)

      real*8, allocatable :: co(:,:),rV(:),cpoints(:,:)

      real*8, allocatable :: dis(:,:)

      character*2, allocatable :: lab(:)

      integer, allocatable :: iZ(:),neigh(:,:),nc(:),ksign(:)

      integer, allocatable :: natom1(:)
      real*8, allocatable :: c1(:,:,:),ch1(:,:),d1(:,:,:),q1(:,:,:,:)
      real*8, allocatable :: o1(:,:,:,:,:),h1(:,:,:,:,:,:)

      dimension vanderwaals_radii(110)

      dimension codot(3)
      dimension t1(3),t2(3,3),t3(3,3,3),t4(3,3,3,3)

      character*20 ca,ca1

      integer alab2anum

      data vanderwaals_radii/
     &2.2676727, 2.6456187, 3.4393037, 3.7794547, 3.7794547, 3.2125367,
     &2.9290777, 2.8723857, 2.7778997, 2.9101797, 4.2896807, 3.2692277,
     &3.7794547, 3.9684267, 3.4015087, 3.4015087, 3.3070227, 3.5526877,
     &5.1967497, 3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.7794547,
     &3.7794547, 3.7794547, 3.7794547, 3.0802557, 2.6456187, 2.6267207,
     &3.5337897, 3.7794547, 3.4959957, 3.5904817, 3.4959957, 3.8172487,
     &3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.7794547,
     &3.7794547, 3.7794547, 3.7794547, 3.0802557, 3.2503307, 2.9857687,
     &3.6471737, 4.1007077, 3.7794547, 3.8928377, 3.7416597, 4.0818107,
     &3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.7794547,
     &3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.7794547,
     &3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.7794547,
     &3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.2503307,
     &3.1369477, 2.9290777, 3.7038657, 3.8172487, 3.7794547, 3.7794547,
     &3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.7794547,
     &3.7794547, 3.5148927, 3.7794547, 3.7794547, 3.7794547, 3.7794547,
     &3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.7794547,
     &3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.7794547,
     &3.7794547, 3.7794547/

      bohr=1.d0/1.8897259886d0

      rP=1.4d0
      rho=1.d0
      open(unit=1,file='IN_SOLSURF',status='old')
      read(1,*)rP
      read(1,*)rho
      close(unit=1)


      open(unit=1,file='name.xyz',status='old')
      read(1,*)nat
      allocate(co(nat,3))
      allocate(lab(nat))
      allocate(iZ(nat))
      allocate(rV(nat))
      allocate(nc(nat))
      allocate(neigh(nat,500))
      allocate(dis(nat,nat))

      read(1,*)
      do n=1,nat
       read(1,*)lab(n),(co(n,k),k=1,3)
      enddo
      close(unit=1)

      do n=1,nat
       iZ(n)=alab2anum(lab(n))
       rV(n)=vanderwaals_radii(Iz(n))*bohr
      enddo


c find the "close" neighbours of each atom
      do n=1,nat-1
      do m=n+1,nat
       dis(n,m)=0.d0
       do k=1,3
        dis(n,m)=dis(n,m)+(co(n,k)-co(m,k))**2
       enddo
       dis(n,m)=sqrt(dis(n,m))
       dis(m,n)=dis(n,m)
      enddo
      enddo

      do n=1,nat
      rlim=(rV(n)+rP)*2.d0
      ncount=0
      do m=1,nat
       if(m.eq.n)go to 1
       if(dis(n,m).lt.rlim)then
        ncount=ncount+1
        if(ncount.gt.500)then
        write(6,*)' ncount too big'
        stop
        endif
        neigh(n,ncount)=m
       endif
1     enddo
      nc(n)=ncount
 
      enddo

      pi=2.d0*acos(0.d0)
      ndim=0.d0
      do n=1,nat
       rs=rV(n)+rP
       sa=4.d0*pi*rs**2
       npts=nint(sa*rho)
       if(npts.gt.ndim)ndim=npts
      enddo
      allocate(cpoints(ndim,3))

c call the subroutine that puts points on a sphere
      ndot=0
      do n=1,nat
       rs=rV(n)+rP
       call spherepoint(ndim,rs,rho,npts,cpoints)
       do i=1,npts
       do k=1,3
        cpoints(i,k)=cpoints(i,k)+co(n,k)
       enddo
       enddo
c only add points not closer to another atom
       rs2=rs*rs
       do j=1,npts
        do m=1,nc(n)
         m1=neigh(n,m)
         ds1=0.d0
         do k=1,3
          ds1=ds1+(cpoints(j,k)-co(m1,k))**2
         enddo
         if(ds1.lt.rs2)go to 2
        enddo
        ndot=ndot+1

        write(2,101)(cpoints(j,k),k=1,3) 
c        write(2,100)'H ',(codot(k),k=1,3)
2      continue
c end the j loop
       enddo
c end the n loop
      enddo   
100   format(a2,3f13.6)
101   format(3f14.6,e15.6)
      close(unit=2)
      open(unit=2,file='fort.2',status='old')
      open(unit=1,file='SOLACCSURFACE',status='unknown')
      write(1,*)nat+ndot
      write(1,*)' molecule coordinates + solvent acc surface'
      do n=1,nat
       write(1,100)lab(n),(co(n,k),k=1,3)
      enddo
      do n=1,ndot
       read(2,*)(codot(k),k=1,3)
       write(1,100)'He',(codot(k),k=1,3)
      enddo
      close(unit=1)
      close(unit=2)


      write(6,*)' The solvent accessible surface has been evaluated.'
      write(6,*)' The surface is defined by the Van der Waals radii'
      write(6,*)' of the atoms in the molecule, and a given radius'
      write(6,*)' of the solvent molecule.'
      write(6,*)' The solvent molecule radius has been input as '
      write(6,102) rp,' Angstrom'
102   format(f6.2,a9)
      write(6,*)' The surface is described by a set of',ndot,' points,'
      write(6,*)' based on a requested points density'
      write(6,103)' (per square Angstrom) of ',rho
103   format(a26,f10.2)
      write(6,*)' The coordinates are in Angstrom '
      write(6,*)' These points are contained in the file SOLACCSURFACE'
      write(6,*)
104   format(2e15.6)


      end

      subroutine spherepoint(ndim,rs,rho,npts,cpoints)
      implicit double precision(a-h,o-z)
      dimension cpoints(ndim,3)
      pi=2.d0*acos(0.d0)
      sa=4.d0*pi*rs**2
      npts=nint(sa*rho)

      s=3.6d0/sqrt(dfloat(npts))
      rlong=0.d0
      dz=2.d0/dfloat(npts)
      z=1.d0-0.5d0*dz
      do j=1,npts
       r=sqrt(1.d0-z*z)
       cpoints(j,1)=r*cos(rlong)
       cpoints(j,2)=r*sin(rlong)
       cpoints(j,3)=z
       z=z-dz
       rlong=rlong+s/r
      enddo
      do j=1,npts
      do k=1,3
       cpoints(j,k)=cpoints(j,k)*rs
      enddo
      enddo

      return
      end


      integer function alab2anum(alab)

c given the atom label this returns the atomic number

      character*2 :: alab,atomic_labels(110)

      data atomic_labels/'H','He','Li','Be','B','C','N','O','F','Ne',
     .'Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr',
     . 'Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb',
     . 'Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',
     . 'Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu',
     . 'Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os',
     . 'Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac',
     .   'Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No'
     . ,'Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds'/

      integer n

      do n=1,110
        if(atomic_labels(n).eq.alab)then
      alab2anum=n
      return
        endif
      enddo
      alab2anum=0

      return
      end


      subroutine Tfactorfast4(x1,x2,t0,t1,t2,t3,t4)
      implicit double precision(a-h,o-z)

      dimension x1(3),x2(3),r(3),dx(3)
      dimension t1(3),t2(3,3),t3(3,3,3),t4(3,3,3,3)

      dimension rr02(3),rr03(3),rr05(3),rr07(3),rr09(3),rr011(3)
      dimension rr023(3),rr025(3),rr027(3),rr029(3),rr0211(3)

c x1 and x2 are the coordinates of the two atoms

c r is the vector from x2 to x1
c and r0 is the distance, |r|.

      r0=0.d0
      do k=1,3
      r(k)=x1(k)-x2(k)
      r0=r0+r(k)**2
      enddo

      r02=1.d0/r0
      r0=1.d0/sqrt(r0)
      r03=r0*r02
      r05=3.d0*r03*r02
      r07=-5.d0*r05*r02

c  First we equate the t factors, t1, t2, etc to the part ofeach
c that has no delta functions

      do k1=1,3
       rr02(k1)=r(k1)*r02
       rr023(k1)=-3.d0*rr02(k1)
       rr025(k1)=-5.d0*rr02(k1)
       rr027(k1)=-7.d0*rr02(k1)
       rr029(k1)=-9.d0*rr02(k1)
       rr0211(k1)=-11.d0*rr02(k1)
       rr03(k1)=-r(k1)*r03
       rr05(k1)=-3.d0*rr03(k1)*r02
       rr07(k1)=-5.d0*rr05(k1)*r02
       rr09(k1)=-7.d0*rr07(k1)*r02
       rr011(k1)=-9.d0*rr09(k1)*r02
      enddo

      t0=r0

      do k1=1,3
      t1(k1)=rr03(k1)
      do k2=k1,3
      t2(k1,k2)=t1(k1)*rr023(k2)
      do k3=k2,3
      t3(k1,k2,k3)=t2(k1,k2)*rr025(k3)
      do k4=k3,3
      t4(k1,k2,k3,k4)=t3(k1,k2,k3)*rr027(k4)
      enddo
      enddo
      enddo
      enddo

c Now we add on the parts that have one delta function, by using the
c same index for two array indices.

      do k1=1,3

      t2(k1,k1)=t2(k1,k1)-r03
      t3(k1,k1,k1)=t3(k1,k1,k1)+rr05(k1)
c 14
      t4(k1,k1,k1,k1)=t4(k1,k1,k1,k1)+rr07(k1)*r(k1) 

      do k2=k1,3
      t3(k1,k1,k2)=t3(k1,k1,k2)+rr05(k2)
      t3(k1,k2,k2)=t3(k1,k2,k2)+rr05(k1)
c 13
      t4(k1,k1,k1,k2)=t4(k1,k1,k1,k2)+rr07(k1)*r(k2)
c 24
      t4(k1,k2,k2,k2)=t4(k1,k2,k2,k2)+rr07(k1)*r(k2)

      do k3=k2,3
c 12
      t4(k1,k1,k2,k3)=t4(k1,k1,k2,k3)+rr07(k2)*r(k3)
c 23
      t4(k1,k2,k2,k3)=t4(k1,k2,k2,k3)+rr07(k1)*r(k3)
c 34
      t4(k1,k2,k3,k3)=t4(k1,k2,k3,k3)+rr07(k1)*r(k2)

      enddo
      enddo
      enddo


3002  continue

c  Now we add on the parts that have two delta functions in a similar
c  way.

      do k1=1,3
c 13 24 + 14 23
      t4(k1,k1,k1,k1)=t4(k1,k1,k1,k1)+2.d0*r05
      do k2=k1,3
c 12 34
      t4(k1,k1,k2,k2)=t4(k1,k1,k2,k2)+r05
      enddo
      enddo

      return
      end
      subroutine ct1d(c1,d2,t1,ect1d)

      implicit double precision(a-h,o-z)

      dimension d2(3),t1(3)

c  this subroutine evaluates the quadrupole-quadrupole
c  electrostatic interaction

      ect1d=c1*(t1(1)*d2(1)+t1(2)*d2(2)+t1(3)*d2(3))

      return
      end

      subroutine ct2q(c1,q2,t2,ect2q)

      implicit double precision(a-h,o-z)

      dimension q2(3,3),t2(3,3)

c  this subroutine evaluates the quadrupole-quadrupole
c  electrostatic interaction

      dimension a(6),b(10)

      a(1)=t2(1,1)
      a(2)=t2(1,2)
      a(3)=t2(1,3)
      a(4)=t2(2,2)
      a(5)=t2(2,3)
      a(6)=t2(3,3)

      b(1)=q2(1,1)
      b(2)=2.d0*q2(1,2)
      b(3)=2.d0*q2(1,3)
      b(4)=q2(2,2)
      b(5)=2.d0*q2(2,3)
      b(6)=q2(3,3)

      ect2q=c1*(a(1)*b(1)+a(2)*b(2)+a(3)*b(3)+a(4)*b(4)+
     .      a(5)*b(5)+a(6)*b(6))

      return
      end

      subroutine ct3o(c1,o2,t3,ect3o)

      implicit double precision(a-h,o-z)

      dimension o2(3,3,3),t3(3,3,3)

c  this subroutine evaluates the quadrupole-quadrupole
c  electrostatic interaction

      dimension a(10),b(10)

      a(1)=o2(1,1,1)*t3(1,1,1)
      a(2)=o2(1,1,2)*t3(1,1,2)*3.d0
      a(3)=o2(1,1,3)*t3(1,1,3)*3.d0
      a(4)=o2(1,2,2)*t3(1,2,2)*3.d0
      a(5)=o2(1,2,3)*t3(1,2,3)*6.d0
      a(6)=o2(1,3,3)*t3(1,3,3)*3.d0
      a(7)=o2(2,2,2)*t3(2,2,2)
      a(8)=o2(2,2,3)*t3(2,2,3)*3.d0
      a(9)=o2(2,3,3)*t3(2,3,3)*3.d0
      a(10)=o2(3,3,3)*t3(3,3,3)

      ect3o=c1*(a(1)+a(2)+a(3)+a(4)+a(5)+a(6)+
     .          a(7)+a(8)+a(9)+a(10))

      return
      end

      subroutine ct4h(c1,h2,t4,ect4h)

      implicit double precision(a-h,o-z)

      dimension h2(3,3,3,3),t4(3,3,3,3)

c  this subroutine evaluates the quadrupole-quadrupole
c  electrostatic interaction

      dimension a(15),b(15)

      a(1)=t4(1,1,1,1)*h2(1,1,1,1)
      a(2)=t4(1,1,1,2)*h2(1,1,1,2)*4.d0
      a(3)=t4(1,1,1,3)*h2(1,1,1,3)*4.d0
      a(4)=t4(1,1,2,2)*h2(1,1,2,2)*6.d0
      a(5)=t4(1,1,2,3)*h2(1,1,2,3)*12.d0
      a(6)=t4(1,1,3,3)*h2(1,1,3,3)*6.d0
      a(7)=t4(1,2,2,2)*h2(1,2,2,2)*4.d0
      a(8)=t4(1,2,2,3)*h2(1,2,2,3)*12.d0
      a(9)=t4(1,2,3,3)*h2(1,2,3,3)*12.d0
      a(10)=t4(1,3,3,3)*h2(1,3,3,3)*4.d0
      a(11)=t4(2,2,2,2)*h2(2,2,2,2)
      a(12)=t4(2,2,2,3)*h2(2,2,2,3)*4.d0
      a(13)=t4(2,2,3,3)*h2(2,2,3,3)*6.d0
      a(14)=t4(2,3,3,3)*h2(2,3,3,3)*4.d0
      a(15)=t4(3,3,3,3)*h2(3,3,3,3)

      ect4h=c1*(a(1)+a(2)+a(3)+a(4)+a(5)+a(6)+
     .          a(7)+a(8)+a(9)+a(10)+a(11)+a(12)+
     .          a(13)+a(14)+a(15))


      return
      end

      subroutine filelabel(k,ca)
      implicit double precision(a-h,o-z)

c this subroutines returns a character variable ca
c for the integer k

      character*20 ca
      character*1 I0(0:10),anum(10),ca1,ca2,ca3,ca4,ca5

      I0(10)='0'
      I0(0)='0'
      I0(1)='1'
      I0(2)='2'
      I0(3)='3'
      I0(4)='4'
      I0(5)='5'
      I0(6)='6'
      I0(7)='7'
      I0(8)='8'
      I0(9)='9'

      if(k.le.9)then
       ca=I0(k)
      endif

      if(k.gt.9.and.k.le.99)then
      k1=k/10
      k2=k-k1*10
      if(k2.eq.0)k2=10
      ca1=I0(k1)
      ca2=I0(k2)
      ca=ca1//ca2
      endif

      if(k.gt.99.and.k.le.999)then
      k1=k/100
      k3=k/10 - INT(k1)*10
      k4=k - INT(k1)*100 - INT(k3)*10
      ca1=I0(k1)
      ca2=I0(k3)
      ca3=I0(k4)
      ca=ca1//ca2//ca3
      endif

      if(k.gt.999.and.k.le.9999)then
      k1=k/1000
      k3=k/100 - INT(k1)*10
      k4=k/10 - INT(k1)*100 - INT(k3)*10
      k5=k-INT(k1)*1000-k3*100-k4*10
      ca1=I0(k1)
      ca2=I0(k3)
      ca3=I0(k4)
      ca4=I0(k5)
      ca=ca1//ca2//ca3//ca4
      endif

      if(k.gt.9999.and.k.le.99999)then
      k1=k/10000
      k3=k/1000 - INT(k1)*10
      k4=k/100 - INT(k1)*100 - INT(k3)*10
      k5=k/10-INT(k1)*1000-k3*100-k4*10
      k6=k-INT(k1)*10000-k3*1000-k4*100-k5*10
      ca1=I0(k1)
      ca2=I0(k3)
      ca3=I0(k4)
      ca4=I0(k5)
      ca5=I0(k6)
      ca=ca1//ca2//ca3//ca4//ca5
      endif

      return
      end

      subroutine getfragdata(nmax)
      implicit double precision(a-h,o-z)

      integer, allocatable :: nat(:)

      open(unit=3,file='frags.out_Lev1',status='old')
      read(3,*)
      read(3,*) nf
      allocate(nat(nf))
      read(3,*)
      read(3,*)
      read(3,*)
      do n=1,nf
       read(3,*)nat(n)
      enddo
      read(3,*)
      do n=1,nf
       read(3,*)(i1,i=1,nat(n))
      enddo
      read(3,*)
1     continue
      read(3,*,end=2)m
      nat(m)=nat(m)+1
      go to 1
2     continue

      nmax=0
      do n=1,nf
       if(nat(n).gt.nmax)nmax=nat(n)
      enddo

      return
      end

      subroutine symm(qu,oct,hex)
      implicit double precision(a-h,o-z)

      dimension qu(3,3),oct(3,3,3),hex(3,3,3,3)

c symmetrize

      do i=1,3
      do j=i,3
       qu(j,i)=qu(i,j)
      enddo
      enddo

      do i=1,3
      do j=i,3
      do k=j,3
      oct(j,i,k)=oct(i,j,k)
      oct(k,j,i)=oct(i,j,k)
      oct(i,k,j)=oct(i,j,k)
      oct(k,i,j)=oct(i,j,k)
      oct(j,k,i)=oct(i,j,k)
      enddo
      enddo
      enddo

      do k1=1,3
      do k2=k1,3
      do k3=k2,3
      do k4=k3,3
      hex(k2,k1,k3,k4)=hex(k1,k2,k3,k4)
      hex(k3,k2,k1,k4)=hex(k1,k2,k3,k4)
      hex(k4,k2,k3,k1)=hex(k1,k2,k3,k4)
      hex(k1,k3,k2,k4)=hex(k1,k2,k3,k4)
      hex(k1,k4,k3,k2)=hex(k1,k2,k3,k4)
      hex(k1,k2,k4,k3)=hex(k1,k2,k3,k4)

      hex(k2,k1,k4,k3)=hex(k1,k2,k3,k4)
      hex(k3,k4,k1,k2)=hex(k1,k2,k3,k4)
      hex(k4,k3,k2,k1)=hex(k1,k2,k3,k4)

      hex(k1,k4,k2,k3)=hex(k1,k2,k3,k4)
      hex(k1,k3,k4,k2)=hex(k1,k2,k3,k4)

      hex(k3,k2,k4,k1)=hex(k1,k2,k3,k4)
      hex(k4,k2,k1,k3)=hex(k1,k2,k3,k4)

      hex(k2,k4,k3,k1)=hex(k1,k2,k3,k4)
      hex(k4,k1,k3,k2)=hex(k1,k2,k3,k4)

      hex(k2,k3,k1,k4)=hex(k1,k2,k3,k4)
      hex(k3,k1,k2,k4)=hex(k1,k2,k3,k4)

      hex(k2,k3,k4,k1)=hex(k1,k2,k3,k4)
      hex(k4,k1,k2,k3)=hex(k1,k2,k3,k4)

      hex(k3,k4,k2,k1)=hex(k1,k2,k3,k4)

      hex(k2,k4,k1,k3)=hex(k1,k2,k3,k4)
      hex(k3,k1,k4,k2)=hex(k1,k2,k3,k4)
      hex(k4,k3,k1,k2)=hex(k1,k2,k3,k4)

      enddo
      enddo
      enddo
      enddo

      return
      end
