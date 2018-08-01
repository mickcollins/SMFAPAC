      subroutine cancelL1L1(ndim,nf,matstore,nsign)
      implicit double precision(a-h,o-z)

      dimension nsign(ndim),matstore(ndim,2)

      write(6,*)'Entered cancelL1L1'
c order the L1L1 fragments in each interaction
      do k=1,nf
       call piksrt(2,matstore(k,:))
      enddo

      do n1=1,nf-1
      if(nsign(n1).eq.0)go to 1
      do n2=n1+1,nf
      if(nsign(n2).eq.0)go to 2
      do j=1,2
       if(matstore(n1,j).ne.matstore(n2,j))go to 2
      enddo
c we have a match
      nsign(n1)=nsign(n1)+nsign(n2)
      nsign(n2)=0
      if(nsign(n1).eq.0)go to 1
2     continue
      enddo
1     continue
      enddo

c elimiate zeros

      ic=0
      do n=1,nf
       if(nsign(n).ne.0)then

        ic=ic+1
        do j=1,2
         matstore(ic,j)=matstore(n,j)
        enddo
        nsign(ic)=nsign(n)
        
       endif
      enddo

      nf=ic

      if(nf.ge.ndim)then
        write(6,*)' cannot get nf below limit'
        write(6,*)' nf = ',nf
        stop
      endif
      write(6,*)' Finished cancelL1L1'
      return
      end
