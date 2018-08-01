      subroutine Tfactorfast(x1,x2,t0,t1,t2,t3,t4,t5,t6)
      implicit double precision(a-h,o-z)

      dimension x1(3),x2(3),r(3),dx(3)
      dimension t1(3),t2(3,3),t3(3,3,3),t4(3,3,3,3)
     .          ,t5(3,3,3,3,3),t6(3,3,3,3,3,3)

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
      r05=3.d0*r03*r02
      r07=-5.d0*r05*r02

      do k1=1,3
      do k2=1,3
       pair(k1,k2)=r(k1)*r(k2)
      enddo
      enddo
 
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
      do k5=k4,3
      t5(k1,k2,k3,k4,k5)=t4(k1,k2,k3,k4)*rr029(k5)
      do k6=k5,3
      t6(k1,k2,k3,k4,k5,k6)=t5(k1,k2,k3,k4,k5)*rr0211(k6)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo

c Now we add on the parts that have one delta function, by using the
c same index for two array indices.

       do k1=1,3
       t5(k1,k1,k1,k1,k1)=t5(k1,k1,k1,k1,k1)+rr09(k1)*pair(k1,k1)
       do k2=k1,3
       t5(k1,k1,k1,k1,k2)=t5(k1,k1,k1,k1,k2)+rr09(k1)*pair(k1,k2)
       t5(k1,k2,k2,k2,k2)=t5(k1,k2,k2,k2,k2)+rr09(k1)*pair(k2,k2)
       do k3=k2,3
       t5(k1,k2,k2,k2,k3)=t5(k1,k2,k2,k2,k3)+rr09(k1)*pair(k2,k3)
       t5(k1,k2,k3,k3,k3)=t5(k1,k2,k3,k3,k3)+rr09(k1)*pair(k2,k3)
       t5(k1,k1,k1,k2,k3)=t5(k1,k1,k1,k2,k3)+rr09(k1)*pair(k2,k3)
       do k4=k3,3
       t5(k1,k2,k3,k3,k4)=t5(k1,k2,k3,k3,k4)+rr09(k1)*pair(k2,k4)
       t5(k1,k2,k2,k3,k4)=t5(k1,k2,k2,k3,k4)+rr09(k1)*pair(k3,k4)
       t5(k1,k1,k2,k3,k4)=t5(k1,k1,k2,k3,k4)+rr09(k2)*pair(k3,k4)
       t5(k1,k2,k3,k4,k4)=t5(k1,k2,k3,k4,k4)+rr09(k1)*pair(k2,k3)
       enddo
       enddo
       enddo
       enddo


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

      do k1=1,3
      do k2=k1,3
      do k3=k2,3
c 12  34
      t5(k1,k1,k2,k2,k3)=
     .t5(k1,k1,k2,k2,k3)+rr07(k3)
c 23  45
      t5(k1,k2,k2,k3,k3)=
     .t5(k1,k2,k2,k3,k3)+rr07(k1)
c 12  45
      t5(k1,k1,k2,k3,k3)=
     .t5(k1,k1,k2,k3,k3)+rr07(k2)

      enddo
      enddo
      enddo

      do k1=1,3
c 13 25 + 14 25 + 14 35 + 15 23 + 15 24 + 15 34
      t5(k1,k1,k1,k1,k1)=
     .t5(k1,k1,k1,k1,k1)+6.d0*rr07(k1)
      do k2=k1,3
c 13  45
      t5(k1,k1,k1,k2,k2)=
     .t5(k1,k1,k1,k2,k2)+rr07(k1)
c 12  35
      t5(k1,k1,k2,k2,k2)=
     .t5(k1,k1,k2,k2,k2)+rr07(k2)
c 14  23
      t5(k1,k1,k1,k1,k2)=
     .t5(k1,k1,k1,k1,k2)+rr07(k2)
c 24  35
      t5(k1,k2,k2,k2,k2)=
     .t5(k1,k2,k2,k2,k2)+rr07(k1)
c 13  24
      t5(k1,k1,k1,k1,k2)=
     .t5(k1,k1,k1,k1,k2)+rr07(k2)
c 25  34
      t5(k1,k2,k2,k2,k2)=
     .t5(k1,k2,k2,k2,k2)+rr07(k1)


      enddo
      enddo


3001  continue

c  now t6 terms

      go to 3003

c  t6 has a part with three delta functions
c  t6 is left out of this code

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
      write(6,*)'t3'
      do k1=1,3
      do k2=k1,3
      do k3=k2,3
       write(6,3)k1,k2,k3,t3(k1,k2,k3)
      enddo
      enddo
      enddo
      write(6,*)'t4'
      do k1=1,3
      do k2=k1,3
      do k3=k2,3
      do k4=k3,3
       write(6,4)k1,k2,k3,k4,t4(k1,k2,k3,k4)
      enddo
      enddo
      enddo
      enddo
      write(6,*)'t5'
      do k1=1,3
      do k2=k1,3
      do k3=k2,3
      do k4=k3,3
      do k5=k4,3
       write(6,5)k1,k2,k3,k4,k5,t5(k1,k2,k3,k4,k5)
      enddo
      enddo
      enddo
      enddo
      enddo
      write(6,*)'t6'
      do k1=1,3
      do k2=k1,3
      do k3=k2,3
      do k4=k3,3
      do k5=k4,3
      do k6=k5,3
       write(6,6)k1,k2,k3,k4,k5,k6,t6(k1,k2,k3,k4,k5,k6)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo

1     format(i3,e15.6)
2     format(2i3,e15.6)
3     format(3i3,e15.6)
4     format(4i3,e15.6)
5     format(5i3,e15.6)
6     format(6i3,e15.6)

3003  continue
      return
      end





	
