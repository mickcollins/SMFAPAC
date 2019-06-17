      subroutine symm(qu,oct,hex)
      implicit double precision(a-h,o-z)

      dimension qu(3,3),oct(3,3,3),hex(3,3,3,3)

c symmetrize

      do i=1,3
      do j=i,3
       qu(j,i)=qu(i,j)
      enddo
      enddo

      return

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
