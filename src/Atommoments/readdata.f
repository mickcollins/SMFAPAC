      subroutine readdata
      use momentheader

      real*8, allocatable :: ch1(:)
      real*8, allocatable :: d1(:,:),q1(:,:,:),o1(:,:,:,:),h1(:,:,:,:,:)

      character*20 ca,ca1

      bohr=1.d0/1.8897259886d0

      allocate(ch1(natom))
      allocate(d1(natom,3))
      allocate(q1(natom,3,3))
      allocate(o1(natom,3,3,3))
      allocate(h1(natom,3,3,3,3))

c read in the electrostatic data for nL1 frags
      do n=1,nL1
       call filelabel(n,ca)
       n1=index(ca,' ')-1
       ca1='nb.'//ca(1:n1)//'.0.cart'
       open(unit=1,file=ca1,status='old')
       read(1,*)
       read(1,*)natom1
       do i=1,natom1
        read(1,*)
        do k=1,3
         read(1,*)c1
c convert to Bohr
        enddo
        read(1,*)
        read(1,*)ch1(i)
        read(1,*)
        do k=1,3
         read(1,*)d1(i,k)
        enddo
        read(1,*)
        do k1=1,3
        do k2=k1,3
        read(1,*)q1(i,k1,k2)
        enddo
        enddo
        read(1,*)
        do k1=1,3
        do k2=k1,3
        do k3=k2,3
        read(1,*)o1(i,k1,k2,k3)
        enddo
        enddo
        enddo
        read(1,*)
        do k1=1,3
        do k2=k1,3
        do k3=k2,3
        do k4=k3,3
        read(1,*)h1(i,k1,k2,k3,k4)
        enddo
        enddo
        enddo
        enddo
c        call symm(q1(i,:,:),o1(i,:,:,:),h1(i,:,:,:,:))
c end the i loop
       enddo
       close(unit=1)
c assign to real atoms
       do i=1,nat1(n)
       do j=1,num1(n,i)
       ch(natnum1(n,i,j))=ch(natnum1(n,i,j))+isign(n)*ch1(i)*w1(n,i,j)
       do k1=1,3
       d(natnum1(n,i,j),k1)=d(natnum1(n,i,j),k1)+
     .      isign(n)*d1(i,k1)*w1(n,i,j)
       do k2=k1,3
       q(natnum1(n,i,j),k1,k2)=q(natnum1(n,i,j),k1,k2)+
     .      isign(n)*q1(i,k1,k2)*w1(n,i,j)
       do k3=k2,3
       o(natnum1(n,i,j),k1,k2,k3)=o(natnum1(n,i,j),k1,k2,k3)+
     .      isign(n)*o1(i,k1,k2,k3)*w1(n,i,j)
       do k4=k3,3
       h(natnum1(n,i,j),k1,k2,k3,k4)=h(natnum1(n,i,j),k1,k2,k3,k4)+
     .      isign(n)*h1(i,k1,k2,k3,k4)*w1(n,i,j)
       enddo
       enddo
       enddo
       enddo
       enddo
       enddo

c end the n loop
      enddo

c     do n=1,natom
c      call symm(q(n,:,:),o(n,:,:,:),h(n,:,:,:,:))
c     enddo



      deallocate(d1)
      deallocate(q1)
      deallocate(o1)
      deallocate(h1)


      return
      end
