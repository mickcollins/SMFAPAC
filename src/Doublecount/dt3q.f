      subroutine dt3q(d1,q2,t3,edt3q)

      implicit double precision(a-h,o-z)

      dimension d1(3),q2(3,3),t3(3,3,3)

c  this subroutine evaluates the quadrupole-quadrupole
c  electrostatic interaction

      dimension a(10),b(10)

      a(1)=t3(1,1,1)
      a(2)=t3(1,1,2)
      a(3)=t3(1,1,3)
      a(4)=t3(1,2,2)
      a(5)=t3(1,2,3)
      a(6)=t3(1,3,3)
      a(7)=t3(2,2,2)
      a(8)=t3(2,2,3)
      a(9)=t3(2,3,3)
      a(10)=t3(3,3,3)

      b(1)=d1(1)*q2(1,1)
      b(2)=2.d0*d1(1)*q2(1,2)+d1(2)*q2(1,1)
      b(3)=2.d0*d1(1)*q2(1,3)+d1(3)*q2(1,1)
      b(4)=2.d0*d1(2)*q2(1,2)+d1(1)*q2(2,2)
      b(5)=2.d0*(d1(1)*q2(2,3)+d1(2)*q2(1,3)+d1(3)*q2(1,2))
      b(6)=2.d0*d1(3)*q2(1,3)+d1(1)*q2(3,3)
      b(7)=d1(2)*q2(2,2)
      b(8)=2.d0*d1(2)*q2(2,3)+d1(3)*q2(2,2)
      b(9)=2.d0*d1(3)*q2(2,3)+d1(2)*q2(3,3)
      b(10)=d1(3)*q2(3,3)

      edt3q=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)+a(4)*b(4)+
     .      a(5)*b(5)+a(6)*b(6)+a(7)*b(7)+a(8)*b(8)+a(9)*b(9)+
     .      a(10)*b(10)

      return
      end

