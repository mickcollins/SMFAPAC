      subroutine zero(natom1,c1,ch1,d1,q1,o1,h1,n0,c0)
      implicit double precision(a-h,o-z)

      dimension c1(natom1,3),ch1(natom1),d1(natom1,3),q1(natom1,3,3)
      dimension o1(natom1,3,3,3),h1(natom1,3,3,3,3)
      dimension c0(60,3)

      dimension mat(natom1)

      do n=1,natom1
       mat(n)=1
      enddo

c we will count the number of atoms deleted
      ntot=0
      do n=1,natom1

       match=0
       do m=1,n0
        dist=0.d0
        do k=1,3
         dist=dist+(c1(n,k)-c0(m,k))**2
        enddo
        if(dist.lt.0.1d0)then
         match=1
         go to 1
        endif
c end the loop over m
       enddo
1      continue
       if(match.eq.1)then
        mat(n)=0
        ntot=ntot+1
        ch1(n)=0.d0
        do k1=1,3
         d1(n,k1)=0.d0
         do k2=1,3
          q1(n,k1,k2)=0.d0
          do k3=1,3
           o1(n,k1,k2,k3)=0.d0
           do k4=1,3
            h1(n,k1,k2,k3,k4)=0.d0
           enddo
          enddo
         enddo
        enddo
       endif
c end the loop over n
      enddo

c correct the charges on the remaining atoms
      if(ntot.lt.natom1)then

c calc the remaining charge
       totch=0.d0
       do n=1,natom1
       totch=totch+ch1(n)
       enddo
c this should be corrected to the nearest integer
       nch=nint(totch)
       delc=(totch-dble(float(nch)))/dble(float(natom1-ntot))
       do n=1,natom1
        if(mat(n).eq.1)ch1(n)=ch1(n)-delc
       enddo
      endif

      return
      end
