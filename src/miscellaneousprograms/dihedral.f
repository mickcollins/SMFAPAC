      program dihedral
      implicit double precision(a-h,o-z)

      real*8, allocatable :: co(:,:)
      character*2, allocatable :: lab(:)
      character*80 filename

      dimension x1(3),x2(3),x3(3),x4(3)
      dimension r12(3),r23(3),r34(3)
      dimension a(3),b(3),c(3)
      dimension acb(3),acc(3),ccb(3)

      pi=2.d0*acos(0.d0)

      open(unit=1,file='xyzFILENAME',status='old')
      read(1,100)filename
      close(unit=1)
100   format(a80)
      open(unit=1,file=trim(filename),status='old')
      read(1,*)natom
      allocate(lab(natom))
      allocate(co(natom,3))
      read(1,*)
      do n=1,natom
       read(1,*)lab(n),(co(n,k),k=1,3)
      enddo

      read(5,*) n1,n2,n3,n4

      do k=1,3
       x1(k)=co(n1,k)
       x2(k)=co(n2,k)
       x3(k)=co(n3,k)
       x4(k)=co(n4,k)
      enddo

      do k=1,3
       r12(k)=x2(k)-x1(k)
       r23(k)=x3(k)-x2(k)
       r34(k)=x4(k)-x3(k)
       c(k)=r23(k)
      enddo
      a(1)=r34(2)*r23(3)-r34(3)*r23(2)
      a(2)=r34(3)*r23(1)-r34(1)*r23(3)
      a(3)=r34(1)*r23(2)-r34(2)*r23(1)
      b(1)=r23(2)*r12(3)-r23(3)*r12(2)
      b(2)=r23(3)*r12(1)-r23(1)*r12(3)
      b(3)=r23(1)*r12(2)-r23(2)*r12(1)
      ccb(1)=c(2)*b(3)-c(3)*b(2)
      ccb(2)=c(3)*b(1)-c(1)*b(3)
      ccb(3)=c(1)*b(2)-c(2)*b(1)

      acc(1)=a(2)*c(3)-a(3)*c(2)
      acc(2)=a(3)*c(1)-a(1)*c(3)
      acc(3)=a(1)*c(2)-a(2)*c(1)

      acb(1)=a(2)*b(3)-a(3)*b(2)
      acb(2)=a(3)*b(1)-a(1)*b(3)
      acb(3)=a(1)*b(2)-a(2)*b(1)

      sum1=0.d0
      sum2=0.d0
      sum3=0.d0
      sum4=0.d0
      sum5=0.d0
      do k=1,3
       sum1=sum1+a(k)*a(k)
       sum2=sum2+b(k)*b(k)
       sum3=sum3+c(k)*c(k)
       sum4=sum4+a(k)*b(k)
       sum5=sum5+c(k)*acb(k)
      enddo
      sum1=sqrt(sum1)
      sum2=sqrt(sum2)
      sum3=sqrt(sum3)
      sumab=sum1*sum2
      cd=sum4/sumab
      sumabc=sumab*sum3
      sd=-sum5/sumabc

      phi=acos(cd)
      if(sd.lt.0.d0)phi=-phi
      phi=phi*180.d0/pi

      write(6,101)n1,n2,n3,n4,phi
101   format(4i10,f13.6)
      end
