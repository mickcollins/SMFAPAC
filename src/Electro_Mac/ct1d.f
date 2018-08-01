      subroutine ct1d(c1,d2,t1,ect1d)

      implicit double precision(a-h,o-z)

      dimension d2(3),t1(3)

c  this subroutine evaluates the quadrupole-quadrupole
c  electrostatic interaction

      ect1d=c1*(t1(1)*d2(1)+t1(2)*d2(2)+t1(3)*d2(3))

      return
      end

