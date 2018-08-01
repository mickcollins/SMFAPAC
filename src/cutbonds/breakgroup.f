      subroutine breakgroup

      use fractdata
      implicit double precision(a-h,o-z)

      integer, allocatable   :: ma1(:,:)
      allocate(ma1(natom,natom))

c find the group on which to base the eleimination of groups

      nafrag=numat(nfrag)

      ma=0

      do n1=1,nafrag
      do n2=1,nafrag
       ma(n1,n2)=iabs(ilink(n1,n2))
      enddo
      enddo

      do n1=1,nafrag
       ma(n1,n1)=1
      enddo


      if(Level.eq.1) go to 1

      do m=1,Level-1

      do n1=1,nafrag
      do n2=1,nafrag
       ma1(n1,n2)=0
      enddo
      enddo

       do n1=1,nafrag
       do n2=1,nafrag
        if(ma(n1,n2).eq.1)then
         do k=1,itp(n2)+1
          ma1(n1,iabs(ibo(n2,k)))=1
          ma1(iabs(ibo(n2,k)),n1)=1
         enddo
        endif
       enddo
       enddo

      do n1=1,nafrag
      do n2=1,nafrag
       if(ma1(n1,n2).eq.1)ma(n1,n2)=1
      enddo
      enddo

       enddo
1      continue

      if(nocapsatall.eq.0.and.Level.eq.2)then
       call three
       call four
      endif

      go to 666

      if(Level.eq.1)go to 666
      if(nocapsatall.eq.1)go to 666

c check for clashing caps
      do n=1,nafrag

10    continue

c look for possible gap between two
       do n1=1,nafrag
       if(ma(n,n1).eq.0)go to 900
       do n2=1,nafrag
       if(n2.eq.n1)go to 901
       if(ma(n,n2).eq.0)go to 901

        do k1=1,itp(n1)+1
         m1=ibo(n1,k1)
         if(iabs(m1).eq.n2)go to 902
         if(ma(n,iabs(m1)).eq.1)go to 902
        do k2=1,itp(n2)+1
         m2=ibo(n2,k2)
         if(iabs(m2).eq.n1)go to 903
         if(ma(n,iabs(m2)).eq.1)go to 903

         if(iabs(m1).eq.iabs(m2))then
c we have a possible "triangle" clash at m1=m2
          if(m1.lt.0.or.m2.lt.0)then
           go to 903
          else
           ma(n,iabs(m1))=1
           go to 10
          endif
         endif

c is m2 a neighbour of m1?
         do k3=1,itp(iabs(m1))+1
          m3=ibo(iabs(m1),k3)
           if(iabs(m3).eq.iabs(m2))then
c we have a possible "rectangle" clash
            if(m1.lt.0.or.m2.lt.0)then
             go to 904
            else
             ma(n,iabs(m1))=1
             ma(n,iabs(m2))=1
             go to 10
            endif
           endif
904      continue
c end the k3 loop
         enddo

c end the k2 loop
903     enddo
c end the k1 loop
902     enddo
c end the n2 loop
901    enddo
c end the n1 loop
900    enddo
c end the n loop
      enddo

666   continue

c find the first group that is not connected to everything else

      match=0
      do n1=1,nafrag

      isum=0
      do n2=1,nafrag
       isum=isum+ma(n1,n2)
      enddo
      if(isum.lt.nafrag)then
       match=1
       n3=n1
       go to 2
      endif

      enddo
2     continue
      if(match.eq.0)then
       write(6,*)' match = 0 for numat = ',nafrag,nfrag
       nstop(nfrag)=1
      deallocate(ma1)     
       return
      endif

c the reference group is n3
      nelim=n3

      deallocate(ma1)

      return
      end

     
