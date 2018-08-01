      program bond
      implicit double precision(a-h,o-z)

      real*8, allocatable :: c(:,:)
      character*2, allocatable :: lab(:)
      character*80 filename

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

      read(5,*) n1,n2
      r=0.d0
      do k=1,3
        r=r+(c(n1,k)-c(n2,k))**2
      enddo
      write(6,101)n1,n2,sqrt(r)
101   format(2i10,f13.6)
      end
