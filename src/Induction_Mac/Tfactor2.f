      subroutine Tfactor2(x1,x2,t0,t1,t2)
      implicit double precision(a-h,o-z)

      dimension x1(3),x2(3),r(3),dx(3)
      dimension t1(3),t2(3,3)

      dimension rr02(3),rr03(3),rr05(3),rr07(3),rr09(3),rr011(3)
      dimension rr023(3),rr025(3),rr027(3),rr029(3),rr0211(3)

      dimension pair(3,3)
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

c  First we equate the t factors, t1, t2, etc to the part ofeach
c that has no delta functions

      do k1=1,3
       rr02(k1)=r(k1)*r02
       rr023(k1)=-3.d0*rr02(k1)
       rr03(k1)=-r(k1)*r03
      enddo

      t0=r0

      do k1=1,3
      t1(k1)=rr03(k1)
      do k2=k1,3
      t2(k1,k2)=t1(k1)*rr023(k2)
      enddo
      enddo

c Now we add on the parts that have one delta function, by using the
c same index for two array indices.


      do k1=1,3
      t2(k1,k1)=t2(k1,k1)-r03
      enddo

      go to 3003

c temp check
      write(6,*)'t0'
      write(6,*)t0
      write(6,*)'t1'
      do k1=1,3
      write(6,1)k1,t1(k1)
      enddo
      write(6,*)'t2'
      do k1=1,3
      do k2=k1,3
       write(6,2)k1,k2,t2(k1,k2)
      enddo
      enddo

1     format(i3,e15.6)
2     format(2i3,e15.6)

3003  continue
      return
      end





