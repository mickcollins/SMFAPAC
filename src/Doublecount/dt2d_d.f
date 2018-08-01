      subroutine dt2d_d(d1,d2,t3,dedt3d)

      implicit double precision(a-h,o-z)

      dimension d1(3),d2(3),t3(3,3,3),dedt3d(3)

c  this subroutine evaluates the derivatives of the dipole-dipole
c  electrostatic interaction

      dimension a(3,6),b(6)

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

      b(1)=d1(1)*d2(1)
      b(2)=d1(1)*d2(2)+d1(2)*d2(1)
      b(3)=d1(1)*d2(3)+d1(3)*d2(1)
      b(4)=d1(2)*d2(2)
      b(5)=d1(2)*d2(3)+d1(3)*d2(2)
      b(6)=d1(3)*d2(3)

      dedt3d(1)=a(1,1)*b(1)+a(1,2)*b(2)+a(1,3)*b(3)+a(1,4)*b(4)+
     .          a(1,5)*b(5)+a(1,6)*b(6)

      dedt3d(2)=a(2,1)*b(1)+a(2,2)*b(2)+a(2,3)*b(3)+a(2,4)*b(4)+
     .          a(2,5)*b(5)+a(2,6)*b(6)

      dedt3d(3)=a(3,1)*b(1)+a(3,2)*b(2)+a(3,3)*b(3)+a(3,4)*b(4)+
     .          a(3,5)*b(5)+a(3,6)*b(6)

      return
      end

