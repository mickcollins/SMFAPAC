      subroutine ct1d_d(c1,d2,t2,dect1d)

      implicit double precision(a-h,o-z)

      dimension d2(3),t2(3,3),dect1d(3)

c  this subroutine evaluates the quadrupole-quadrupole
c  electrostatic interaction

      dect1d(1)=c1*(t2(1,1)*d2(1)+t2(1,2)*d2(2)+t2(1,3)*d2(3))
      dect1d(2)=c1*(t2(1,2)*d2(1)+t2(2,2)*d2(2)+t2(2,3)*d2(3))
      dect1d(3)=c1*(t2(1,3)*d2(1)+t2(2,3)*d2(2)+t2(3,3)*d2(3))

      return
      end

