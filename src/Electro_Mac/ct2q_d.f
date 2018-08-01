      subroutine ct2q_d(c1,q2,t3,dect2q)

      implicit double precision(a-h,o-z)

      dimension q2(3,3),t3(3,3,3),dect2q(3)

c  this subroutine evaluates the quadrupole-quadrupole
c  electrostatic interaction

      dimension a(3,6),b(10)

      a(1,1)=t3(1,1,1)
      a(1,2)=t3(1,1,2)
      a(1,3)=t3(1,1,3)
      a(1,4)=t3(1,2,2)
      a(1,5)=t3(1,2,3)
      a(1,6)=t3(1,3,3)

      a(2,1)=t3(1,1,2)
      a(2,2)=t3(1,2,2)
      a(2,3)=t3(1,2,3)
      a(2,4)=t3(2,2,2)
      a(2,5)=t3(2,2,3)
      a(2,6)=t3(2,3,3)

      a(3,1)=t3(1,1,3)
      a(3,2)=t3(1,2,3)
      a(3,3)=t3(1,3,3)
      a(3,4)=t3(2,2,3)
      a(3,5)=t3(2,3,3)
      a(3,6)=t3(3,3,3)

      b(1)=q2(1,1)
      b(2)=2.d0*q2(1,2)
      b(3)=2.d0*q2(1,3)
      b(4)=q2(2,2)
      b(5)=2.d0*q2(2,3)
      b(6)=q2(3,3)

      do k=1,3
      dect2q(k)=
     .c1*(a(k,1)*b(1)+a(k,2)*b(2)+a(k,3)*b(3)+a(k,4)*b(4)+
     .      a(k,5)*b(5)+a(k,6)*b(6))
      enddo

      return
      end

