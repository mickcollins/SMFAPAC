      program angle
      implicit double precision(a-h,o-z)

      real*8, allocatable :: c(:,:)
      character*2, allocatable :: lab(:)
      character*80 filename

      pi=2.d0*acos(0.d0)

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

      read(5,*) n1,n2,n3
      r1=0.d0
      r2=0.d0
      r3=0.d0
      do k=1,3
        r1=r1+(c(n1,k)-c(n2,k))**2
        r2=r2+(c(n3,k)-c(n2,k))**2
        r3=r3+(c(n1,k)-c(n3,k))**2
      enddo
      anom=r1+r2-r3
      denom=2.d0*sqrt(r1*r2)
      ang=anom/denom
      ang=acos(ang)*180.d0/pi
      write(6,101)n1,n2,n3,ang
101   format(3i10,f13.6)
      end
