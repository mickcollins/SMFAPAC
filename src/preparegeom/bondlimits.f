      subroutine bondlimits(b)
      use AtomicData
      implicit double precision(a-h,o-z)

      dimension b(110,110)

      do n1=1,91
      do n2=n1,91
       b(n1,n2)=arad(n1)+arad(n2)+0.08d0
      enddo
      enddo

      b(7,33)=1.800d0
      b(8,33)=1.690d0
      b(15,33)=2.250d0
      b(16,33)=2.096d0

      b(6,6)=1.440d0
      b(6,7)=1.410d0
      b(6,8)=1.412d0
      b(6,15)=1.757d0
      b(6,16)=1.735d0
      b(6,52)=2.082d0

      b(7,7)=1.400d0
      b(7,8)=1.300d0
      b(7,15)=1.620d0
      b(7,16)=1.572d0
      b(7,34)=1.794d0

      b(8,15)=1.619d0
      b(8,16)=1.530d0
      b(8,34)=1.770d0

      b(15,15)=2.122d0
      b(15,16)=1.98d0
      b(15,34)=2.150d0
      b(15,52)=2.400d0

      do n1=1,91
      do n2=n1,91
       b(n2,n1)=b(n1,n2)
      enddo
      enddo

      return
      end
