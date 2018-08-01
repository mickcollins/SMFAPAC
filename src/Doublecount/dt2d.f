      subroutine dt2d(d1,d2,t2,edt2d)

      implicit double precision(a-h,o-z)

      dimension d1(3),d2(3),t2(3,3)

c  this subroutine evaluates the quadrupole-quadrupole
c  electrostatic interaction

      dimension a(6),b(6)

      a(1)=t2(1,1)
      a(2)=t2(1,2)
      a(3)=t2(1,3)
      a(4)=t2(2,2)
      a(5)=t2(2,3)
      a(6)=t2(3,3)

      b(1)=d1(1)*d2(1)
      b(2)=d1(1)*d2(2)+d1(2)*d2(1)
      b(3)=d1(1)*d2(3)+d1(3)*d2(1)
      b(4)=d1(2)*d2(2)
      b(5)=d1(2)*d2(3)+d1(3)*d2(2)
      b(6)=d1(3)*d2(3)

      edt2d=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)+a(4)*b(4)+
     .      a(5)*b(5)+a(6)*b(6)

      return
      end

