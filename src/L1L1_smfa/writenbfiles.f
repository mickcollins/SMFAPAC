      subroutine writenbfiles(ngroups,nffinal0,nch0,nat0,lab0,
     .           c0,nat1,lab1,c1,nffinal,nffinal1,nbextra,numcharges,numchggrps,
     .           numat,natstore,nch1,sym,numb)

      implicit double precision(a-h,o-z)

      dimension nch0(ngroups),nat0(ngroups),nch1(2*nffinal1),
     .          natstore(2*nffinal1,8),numat(2*nffinal1),
     .          nat1(2*nffinal1),numchggrps(numcharges)

      dimension c0(ngroups,60,3),c1(2*nffinal1,60,3)

      dimension num(2000),nbg(2000,6)
      dimension ptch(2000),we(2000,6)

      character*2 lab0(ngroups,60),lab1(2*nffinal1,60)

      character*20 ca,ca1
      character*27 ca2

      dimension num0(60),num1(60)

c open the file that will allocate atoms ids
c to these bg charges
       open(unit=4,file='NBbgidentities',status='unknown')
       write(4,*)nffinal+nbextra

c write to fort.23 which will be added to OUT_ELECTRONS
       write(23,*)' The number of Level 1 fragments is'
       write(23,*)nffinal+nbextra
       write(23,*)'The numbers of atoms and electrons in each fragment'

      do n=1,nffinal+nbextra
       call filelabel(n,ca)
       n1=index(ca,' ')-1
       ca1='nb.'//ca(1:n1)//'.0.com'
       open(unit=1,file=ca1,status='unknown')
c calculate the charge and multiplicity
       if(numat(n).eq.1)then
        ncharge=nch0(natstore(n,1))
        call Atnum(60,nat0(natstore(n,1)),lab0(natstore(n,1),:),
     .            num0,sym,numb,nflag)
        if(nflag.eq.1)then
c        write(6,*)' n ',n
c        write(6,*)' nat0(natstore(n,1)) ',nat0(natstore(n,1))
         do i=1,nat0(natstore(n,1))
c         write(6,*)i,lab0(natstore(n,1),i)
         enddo
         stop
        endif
        ne=-nch0(natstore(n,1))
        do i=1,nat0(natstore(n,1))
         ne=ne+num0(i)
        enddo
        ne2=ne/2
        ne3=ne-2*ne2
        if(ne3.eq.0)then
         multiplicity=1
        else
         multiplicity=2
        endif

       else
        ncharge=nch1(n)
        call Atnum(60,nat1(n),lab1(n,:),
     .            num1,sym,numb,nflag)
        if(nflag.eq.1)then
c        write(6,*)' n ',n
c        write(6,*)' nat1(n) ',nat1(n)
         do i=1,nat1(n)
c         write(6,*)i,lab1(n,i)
         enddo
         stop
        endif
        ne=-nch1(n)
        do i=1,nat1(n)
         ne=ne+num1(i)
        enddo
        ne2=ne/2
        ne3=ne-2*ne2
        if(ne3.eq.0)then
         multiplicity=1
        else
         multiplicity=2
        endif

       endif


       write(1,*)ncharge,multiplicity
c write the coordinates
       if(numat(n).eq.1)then
        do i=1,nat0(natstore(n,1))
         write(1,100)lab0(natstore(n,1),i),(c0(natstore(n,1),i,k),k=1,3)
        enddo
        netatoms=nat0(natstore(n,1))
       else
        do i=1,nat1(n)
         write(1,100)lab1(n,i),(c1(n,i,k),k=1,3)
        enddo
        netatoms=nat1(n)
       endif
100    format(a2,3f13.6)
c add a blank line
       write(1,*)
c write the number of atoms and electrons for test information
       write(23,*)netatoms,ne,ca1

c skip the embedded charges in this version of L1L1
      go to 3000

c add any point charges if necessary
      if(numcharges.gt.0)then
      jc=0
      do k=1,numcharges
c skip the charge group if it is contained in this L1 frag
        do j=1,numat(n)
         if(natstore(n,j).eq.numchggrps(k))go to 440
        enddo
c this charge group can be added 
       call filelabel(numchggrps(k),ca)
       n1=index(ca,' ')-1
       ca1='charge.'//ca(1:n1)//'.grp'
       ca2='charge.'//ca(1:n1)//'.grp'//'_atoms'
       open(unit=2,file=ca1,status='old')
       open(unit=3,file=ca2,status='old')
441    continue
       read(2,*,end=442)x1,y1,z1,ptch1
       write(1,101)x1,y1,z1,ptch1
       jc=jc+1
       read(3,*)num(jc),(nbg(jc,l),we(jc,l),l=1,num(jc))
       ptch(jc)=ptch1
       go to 441
442    continue
       close(unit=2)
       close(unit=3)
101    format(4f15.7)

c end the loop over charge files
440   continue
      enddo
c add a blank line
       write(1,*)

      endif
c output the bg charges identities
      if(numcharges.eq.0)then
       write(4,*)0
      else
       write(4,*)jc
       if(jc.gt.0)then
        do l=1,jc
         write(4,668)ptch(l),num(l),(nbg(l,k),we(l,k),k=1,num(l))
        enddo
       endif
      endif
668    format(f10.6,i4,6(i6,f10.6))
     
3000  continue 
c end the loop over L1 frags
      close(unit=1)
      enddo
      close(unit=4)
      return
      end
