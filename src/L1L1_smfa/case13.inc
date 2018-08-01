      if(nsum.eq.1.and.numat(k1)*numat(k2).eq.4)then
      do j1=1,numat(k1)
      do j2=1,numat(k2)
       if(mata(j1,j2).eq.1)then
        mg1=natstore(k1,j1)
        mg2=natstore(k2,j2)
      call match_mg(2*nffinal1,mg1,numat,natstore,nffinal,nbextra,match)
       match1=match
      call match_mg(2*nffinal1,mg2,numat,natstore,nffinal,nbextra,match)
       match2=match
c fill in the fragments
cc       if(ic.ge.10*nffinal0)then
cc        write(6,*)' ic too big'
cc      call cancelL1L1(10*nffinal0,ic,matstore,nsign)
cc       endif
       ic=ic+1
c record the fragment numbers and the net sign
       matstore(ic,1)=match1
       matstore(ic,2)=match2
       nsign(ic)=ksign(k1)*ksign(k2)
      if(matstore(ic,1).eq.0.or.matstore(ic,2).eq.0)then
       write(6,*)'case13'
       write(6,*)k1,k2
       write(6,*)j1,j2
       write(6,*)numat(k1),numat(k2)
       stop
      endif
c     do ip=1,ic
c     if(matstore(ip,1).eq.0.or.matstore(ip,2).eq.0)then
c      write(6,*)'case13'
c      write(6,*)ip,ic
c      stop
c     endif
c     enddo
       go to 111

       endif
c end the j1,j2 loops
       enddo
       enddo
c end the nsum if
      endif

