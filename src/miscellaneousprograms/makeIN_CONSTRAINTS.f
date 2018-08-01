      program makeIN_CONSTRAINTS
      implicit double precision(a-h,o-z)

      integer, allocatable  :: mb(:),nb(:),ia(:),ja(:),ka(:)
      integer, allocatable  :: id(:),jd(:),kd(:),ld(:)
      real*8, allocatable   :: r(:), dr(:), a(:), da(:)
      real*8, allocatable   :: d(:),dd(:)

      character*1 I0(0:10),ca1,ca2,ca3,ca4
      character*20 ca

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


      open(unit=1,file='IN_OPTSCAN', status='old')

      read(1,*)
      read(1,*)nscan
      if(nscan.eq.0)go to 1
      read(1,*)
      read(1,*)nbonds
      if(nbonds.eq.0)go to 2
      read(1,*)
      allocate(mb(nbonds))
      allocate(nb(nbonds))
      allocate(r(nbonds))
      allocate(dr(nbonds))
      do n=1,nbonds
       read(1,*)mb(n),nb(n),r(n),dr(n)
      enddo
2     read(1,*)
      read(1,*)nangles
      if(nangles.eq.0)go to 3

      allocate(ia(nangles))
      allocate(ja(nangles))
      allocate(ka(nangles))
      allocate(a(nangles))
      allocate(da(nangles))
      read(1,*)
      do n=1,nangles
       read(1,*)ia(n),ja(n),ka(n),a(n),da(n)
      enddo
3     continue
      read(1,*)
      read(1,*)ndiheds
      if(ndiheds.eq.0)go to 44
      allocate(id(ndiheds))
      allocate(jd(ndiheds))
      allocate(kd(ndiheds))
      allocate(ld(ndiheds))
      allocate(d(ndiheds))
      allocate(dd(ndiheds))
      read(1,*)
      do n=1,ndiheds
       read(1,*)id(n),jd(n),kd(n),ld(n),d(n),dd(n)
      enddo
44     continue

      close(unit=1)

      do k=1,nscan+1

c generate the group file name

      if(k.le.9)then
       ca='IN_CONSTRAINTS_'//I0(k)
      endif
      if(k.gt.9.and.k.le.99)then
      k1=k/10
      k2=k-k1*10
      ca1=I0(k1)
      ca2=I0(k2)
      ca='IN_CONSTRAINTS_'//ca1//ca2

      endif
      if(k.gt.99.and.k.le.999)then
      k1=k/100
      k3=k/10 - INT(k1)*10
      k4=k - INT(k1)*100 - INT(k3)*10
      ca1=I0(k1)
      ca2=I0(k3)
      ca3=I0(k4)
      ca='IN_CONSTRAINTS_'//ca1//ca2//ca3
      endif
      if(k.gt.999.and.k.le.9999)then
      k1=k/1000
      k2=k/100 - INT(k1)*10
      k3=k/10 - INT(k2)*10-int(k1)*100
      k4=k - INT(k1)*1000 - INT(k3)*10-int(k2)*100
      ca1=I0(k1)
      ca2=I0(k2)
      ca3=I0(k3)
      ca4=I0(k4)
      ca='IN_CONSTRAINTS_'//ca1//ca2//ca3//ca4
      endif
      open(unit=2,file=ca,status='unknown')

c write the data
      write(2,*)' Enter the number of bond constraints'
      write(2,*)nbonds
      if(nbonds.eq.0)go to 4
      write(2,*)' Enter the atoms in the bond and the length (Angstrom)'
      do n=1,nbonds
       write(2,*)mb(n),nb(n),r(n)+(k-1)*dr(n)
      enddo
4     write(2,*)' Enter the number of angle constraints'
      write(2,*)nangles
      if(nangles.eq.0)go to 5
      write(2,*)' Enter the atoms in the angle & angle value (degrees)'
      do n=1,nangles
       write(2,*)ia(n),ja(n),ka(n),a(n)+(k-1)*da(n)
      enddo
5     write(2,*)' Enter the number of dihedral angles to be constrained'
      write(2,*)ndiheds
      if(ndiheds.eq.0)go to 45
      write(2,*)'Enter the 4 atom numbers and dihedral angles (degrees)'
      do n=1,ndiheds
       write(2,*) id(n),jd(n),kd(n),ld(n),d(n)+(k-1)*dd(n)
      enddo
45    continue     
      close(unit=2)

c end the scan loop
      enddo
1     continue
      end

