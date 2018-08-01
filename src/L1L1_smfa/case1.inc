      if(nsum.eq.3.and.mata(1,1).eq.0)then
c three interactions to be created

       mg=natstore(k1,2)
      call match_mg(2*nffinal1,mg,numat,natstore,nffinal,nbextra,match)
c save match
       match1=match
c fill in the fragments
cc       if(ic.ge.10*nffinal0)then
cc        write(6,*)' ic too big'
cc     call cancelL1L1(10*nffinal0,ic,matstore,nsign)
cc       endif
       ic=ic+1
c record the fragment numbers and the net sign
       matstore(ic,1)=match
       matstore(ic,2)=k2
       nsign(ic)=ksign(k1)*ksign(k2)

c now 2nd of 3 interactions
       mg=natstore(k2,2)
      call match_mg(2*nffinal1,mg,numat,natstore,nffinal,nbextra,match)
c save match
       match2=match
c fill in the fragments
cc       if(ic.ge.10*nffinal0)then
cc        write(6,*)' ic too big'
cc      call cancelL1L1(10*nffinal0,ic,matstore,nsign)
cc       endif
       ic=ic+1
c record the fragment numbers and the net sign
       matstore(ic,1)=k1
       matstore(ic,2)=match
       nsign(ic)=ksign(k1)*ksign(k2)

c now 3rd of three
cc       if(ic.ge.10*nffinal0)then
cc        write(6,*)' ic too big'
cc      call cancelL1L1(10*nffinal0,ic,matstore,nsign)
cc       endif
       ic=ic+1
c record the fragment numbers and the net sign
       matstore(ic,1)=match1
       matstore(ic,2)=match2
       nsign(ic)=-ksign(k1)*ksign(k2)
c     do ip=1,ic
c     if(matstore(ip,1).eq.0.or.matstore(ip,2).eq.0)write(6,*)'case1'
c     enddo
       go to 111
      endif

