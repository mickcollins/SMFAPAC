      subroutine ct2q(c1,q2,t2,ect2q)

      implicit double precision(a-h,o-z)

      dimension q2(3,3),t2(3,3)

c  this subroutine evaluates the quadrupole-quadrupole
c  electrostatic interaction

      dimension a(6),b(10)

      a(1)=t2(1,1)
      a(2)=t2(1,2)
      a(3)=t2(1,3)
      a(4)=t2(2,2)
      a(5)=t2(2,3)
      a(6)=t2(3,3)

      b(1)=q2(1,1)
      b(2)=2.d0*q2(1,2)
      b(3)=2.d0*q2(1,3)
      b(4)=q2(2,2)
      b(5)=2.d0*q2(2,3)
      b(6)=q2(3,3)

      ect2q=c1*(a(1)*b(1)+a(2)*b(2)+a(3)*b(3)+a(4)*b(4)+
     .      a(5)*b(5)+a(6)*b(6))

      return
      end

