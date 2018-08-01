      subroutine filelabel(k,ca)
      implicit double precision(a-h,o-z)

c this subroutines returns a character variable ca
c for the integer k

      character*20 ca
      character*1 I0(0:10),anum(10),ca1,ca2,ca3,ca4,ca5

      I0(10)='0'
      I0(0)='0'
      I0(1)='1'
      I0(2)='2'
      I0(3)='3'
      I0(4)='4'
      I0(5)='5'
      I0(6)='6'
      I0(7)='7'
      I0(8)='8'
      I0(9)='9'

      if(k.le.9)then
       ca=I0(k)
      endif

      if(k.gt.9.and.k.le.99)then
      k1=k/10
      k2=k-k1*10
      if(k2.eq.0)k2=10
      ca1=I0(k1)
      ca2=I0(k2)
      ca=ca1//ca2
      endif

      if(k.gt.99.and.k.le.999)then
      k1=k/100
      k3=k/10 - INT(k1)*10
      k4=k - INT(k1)*100 - INT(k3)*10
      ca1=I0(k1)
      ca2=I0(k3)
      ca3=I0(k4)
      ca=ca1//ca2//ca3
      endif

      if(k.gt.999.and.k.le.9999)then
      k1=k/1000
      k3=k/100 - INT(k1)*10
      k4=k/10 - INT(k1)*100 - INT(k3)*10
      k5=k-INT(k1)*1000-k3*100-k4*10
      ca1=I0(k1)
      ca2=I0(k3)
      ca3=I0(k4)
      ca4=I0(k5)
      ca=ca1//ca2//ca3//ca4
      endif

      if(k.gt.9999.and.k.le.99999)then
      k1=k/10000
      k3=k/1000 - INT(k1)*10
      k4=k/100 - INT(k1)*100 - INT(k3)*10
      k5=k/10-INT(k1)*1000-k3*100-k4*10
      k6=k-INT(k1)*10000-k3*1000-k4*100-k5*10
      ca1=I0(k1)
      ca2=I0(k3)
      ca3=I0(k4)
      ca4=I0(k5)
      ca5=I0(k6)
      ca=ca1//ca2//ca3//ca4//ca5
      endif

      return
      end

