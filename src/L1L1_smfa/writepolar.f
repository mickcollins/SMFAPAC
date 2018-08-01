      subroutine writepolar(ngroups,nch0,nat0,lab0,c0,sym,numb)

      implicit double precision(a-h,o-z)

      dimension nch0(ngroups),nat0(ngroups)

      dimension c0(ngroups,60,3)

      character*2 lab0(ngroups,60),sym(110)

      character*20 ca,ca1

      dimension num(60)

      dimension numb(110)

      write(24,*)' The number of polarisation calcs is'
      write(24,*)ngroups
      write(24,*)' The number of atoms and electrons in each calc is'
      do n=1,ngroups

c calculate the number of electrons
      call Atnum(60,nat0(n),lab0(n,:),num,sym,numb,nflag)
        if(nflag.eq.1)then
         write(6,*)' n ',n
         write(6,*)' nat0(n) ',nat0(n)
         do i=1,nat0(n)
          write(6,*)i,lab0(n,i)
         enddo
         stop
        endif
      ne=nch0(n)
      do i=1,nat0(n)
       ne=ne+num(i)
      enddo

      ne2=ne/2
      ne3=ne-2*ne2
      if(ne3.eq.0)then
       multiplicity=1
      else
       multiplicity=2
      endif
c make the filename
       call filelabel(n,ca)
       n1=index(ca,' ')-1
       ca1='nb.'//ca(1:n1)//'.0-polar.com'
       open(unit=1,file=ca1,status='unknown')
       write(1,*)nch0(n),multiplicity     
       do i=1,nat0(n)
        write(1,100)lab0(n,i),(c0(n,i,k),k=1,3)
       enddo
100    format(a2,3f13.6)
c add a blank line
       write(1,*)
c end the loop over groups
       close(unit=1)
c write out the number of electrons for information
c     write(6,*)'writing 24'
      write(24,*)nat0(n),ne,ca1

      enddo

      return
      end
