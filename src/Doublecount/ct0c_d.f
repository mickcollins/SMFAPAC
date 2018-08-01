      subroutine ct0c_d(c1,c2,t1,dct0c)

      implicit double precision(a-h,o-z)

      dimension t1(3),dct0c(3)

c  this subroutine evaluates thederivative of the charge-charge
c  electrostatic interaction

      c=c1*c2

      do k=1,3
       dct0c(k)=c*t1(k)
      enddo

      return
      end

