      program convertgamhessian
      implicit double precision(a-h,o-z)

      real*8, allocatable :: f(:,:)

      dimension f1(3),f2(3)

      read(5,*)natom
      allocate(F(3*natom,3*natom))

      natom2=natom/2
      ndiff=natom-2*natom2

      do i=1,natom-1,2
       i1=i+1
       n2=3*(i-1)
       n3=3*(i1-1)
       do j=i,natom
        do k=1,3
         n1=3*(j-1)+k
          read(5,100)(f1(k2),k2=1,3),(f2(k3),k3=1,3)
100      format(20x,6f9.6)
         do k2=1,3
          F(n1,n2+k2)=f1(k2)
         enddo
         do k3=1,3
          F(n1,n3+k3)=f2(k3)
         enddo
        enddo
       enddo
      enddo
      if(ndiff.eq.1)then
      i=natom
      n2=3*(i-1)
      j=natom
        do k=1,3
         n1=3*(j-1)+k
          read(5,100)(f1(k2),k2=1,3)
         do k2=1,3
          F(n1,n2+k2)=f1(k2)
         enddo
        enddo
      endif

      n3=3*natom
      do i=1,n3
      do j=1,i
       write(6,101)F(i,j)
      enddo
      enddo
101   format(1x,f9.6)
      end

