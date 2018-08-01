      subroutine project(fcm,c,nat)
      implicit double precision (a-h,o-z)

C     uses the Miller,Handy and Adams algorithm
C     to project out rotational and translational
C     components of teh force constant matrix
C     (strickly only valid at a stationary point
C      but can be useful in geometry optimisations)

C     fcm = force constant matrix
C     C = cartesian geometry
C     Nat = number of atoms

C     Caution - overwrites original force constant matrix
C     fcm with projected version

      dimension c(3,nat),fcm(nat*3,nat*3)
c     temporary arrays
c     dimension f1(nat*3,nat*3),al(nat*3,nat*3)

      real*8, allocatable :: f1(:,:),al(:,:)

      nat3=nat*3
      allocate(f1(nat3,nat3))
      allocate(al(nat3,nat3))

      call zero(al,nat3*6)
      DO  N=1,NAT
      AL(3*N-2,1)=1.0D0
      AL(3*N-1,2)=1.0D0
      AL(3*N,3)=1.0D0
      AL(3*N-2,4)=-C(2,N)
      AL(3*N-1,4)=C(1,N)
      AL(3*N-1,5)=-C(3,N)
      AL(3*N,5)=C(2,N)
      AL(3*N-2,6)=C(3,N)
      AL(3*N,6)=-C(1,N)
      ENDDO

      call  orthy(al,nat3,6)
      call dgemm('N','T',nat3,nat3,6,-1d0,al,nat3,al,nat3,0d0,f1,nat3)

      do i=1,nat3
      f1(i,i)=f1(i,i)+1d0
      enddo

C
      call dgemm('N','N',nat3,nat3,nat3,1d0,f1,nat3,fcm,nat3,0d0,al,
     1 nat3)
      call dgemm('N','N',nat3,nat3,nat3,1d0,al,nat3,f1,nat3,0d0,fcm,
     1   nat3)

      RETURN
      END
      subroutine orthy(a,m,n)
      implicit double precision (a-h,o-z)
      dimension a(m,n), work(n)

c     orthonormalises the columns of matrix A
c     caution : intended as part of subroutine projec
c     not suitable for general use
c     sets small elements of A to zero
c     small aribitrarily defined as 1d-10

      parameter (small=1d-10)

      do j=1,n
      do i=1,m
      if(dabs(a(i,j)).lt.small)a(i,j)=0d0
      enddo
      enddo

      c=dsqrt(ddot(m,a(1,1),1,a(1,1),1))

      if(c.gt.small)then
       a(1:m,1)=a(1:m,1)/c
      else
       a(1:m,1)=0d0
      endif


      do na=2,n

          do nb=1,na-1
          work(nb)=ddot(m,a(1,na),1,a(1,nb),1)
          enddo

          do nb=1,na-1
          call daxpy(m,-work(nb),a(1,nb),1,a(1,na),1)
          enddo

          do i=1,m
          if(dabs(a(i,na)).lt.small)a(i,na)=0d0
          enddo

          c=dsqrt(ddot(m,a(1,na),1,a(1,na),1))
          if(c.gt.small)then
          a(1:m,na)=a(1:m,na)/c
          else
          a(1:m,na)=0d0
          endif

      enddo
      return
      end
      subroutine zero(a,n)
      double precision a(n)
      do i=1,n
      a(i)=0d0
      enddo
      return
      end

