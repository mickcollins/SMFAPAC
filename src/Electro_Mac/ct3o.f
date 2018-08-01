      subroutine ct3o(c1,o2,t3,ect3o)

      implicit double precision(a-h,o-z)

      dimension o2(3,3,3),t3(3,3,3)

c  this subroutine evaluates the quadrupole-quadrupole
c  electrostatic interaction

      dimension a(10),b(10)

      a(1)=o2(1,1,1)*t3(1,1,1)
      a(2)=o2(1,1,2)*t3(1,1,2)*3.d0
      a(3)=o2(1,1,3)*t3(1,1,3)*3.d0
      a(4)=o2(1,2,2)*t3(1,2,2)*3.d0
      a(5)=o2(1,2,3)*t3(1,2,3)*6.d0
      a(6)=o2(1,3,3)*t3(1,3,3)*3.d0
      a(7)=o2(2,2,2)*t3(2,2,2)
      a(8)=o2(2,2,3)*t3(2,2,3)*3.d0
      a(9)=o2(2,3,3)*t3(2,3,3)*3.d0
      a(10)=o2(3,3,3)*t3(3,3,3)

      ect3o=c1*(a(1)+a(2)+a(3)+a(4)+a(5)+a(6)+
     .          a(7)+a(8)+a(9)+a(10))

      return
      end

