      subroutine writeabfiles(ngroups,nffinal0,nf,nch0,nat0,lab0,
     .           c0,nat1,lab1,c1,nffinal,nffinal1,nbextra,numcharges,numchggrps,
     .           numat,natstore,nch1,sym,numb,k1,k2,intsign)

      implicit double precision(a-h,o-z)

      dimension nch0(ngroups),nat0(ngroups),nch1(2*nffinal1),
     .          natstore(2*nffinal1,8),numat(2*nffinal1),
     .          nat1(2*nffinal1),numchggrps(numcharges)

      dimension c0(ngroups,60,3),c1(2*nffinal1,60,3)

      dimension num(2000),nbg(2000,6)
      dimension ptch(2000),we(2000,6)

      character*2 lab0(ngroups,60),lab1(2*nffinal1,60)

      character*20 ca,ca1,ca2
      character*27 da

      dimension num0(60),num1(60)

      dimension numb(110)
      character*2 sym(110)


c make the file name
       call filelabel(k1,ca1)
       n1=index(ca1,' ')-1
       call filelabel(k2,ca2)
       n2=index(ca2,' ')-1
       ca='ab.'//ca1(1:n1)//'.'//ca2(1:n2)//'.com'
       open(unit=1,file=ca,status='unknown')

       da='ABbgidentities.'//ca1(1:n1)//'.'//ca2(1:n2)
       open(unit=4,file=da,status='unknown')

c calculate the charge and multiplicity
       if(numat(k1).eq.1)then

        ncharge=nch0(natstore(k1,1))
        call Atnum(60,nat0(natstore(k1,1)),lab0(natstore(k1,1),:),
     .            num0,sym,numb,nflag)
        if(nflag.eq.1)then
c        write(6,*)' k1,k2 ',k1,k2
c        write(6,*)' nat0(natstore(k1,1)) ',nat0(natstore(k1,1))
         do i=1,nat0(natstore(k1,1))
c         write(6,*)i,lab0(natstore(k1,1),i)
         enddo
         stop
        endif

        netatoms=nat0(natstore(k1,1))
        ne=-nch0(natstore(k1,1))
        do i=1,nat0(natstore(k1,1))
         ne=ne+num0(i)
        enddo

       else
        ncharge=nch1(k1)
        call Atnum(60,nat1(k1),lab1(k1,:),
     .            num1,sym,numb,nflag)
        if(nflag.eq.1)then
c        write(6,*)' k1,k2 ',k1,k2
c        write(6,*)' nat1(k1) ',nat1(k1)
         do i=1,nat1(k1)
c         write(6,*)i,lab1(k1,i),num1(i)
         enddo
         stop
        endif
        netatoms=nat1(k1)
        ne=-nch1(k1)
        do i=1,nat1(k1)
         ne=ne+num1(i)
        enddo


       endif
c repeat for k2
       if(numat(k2).eq.1)then

        ncharge=ncharge+nch0(natstore(k2,1))
        call Atnum(60,nat0(natstore(k2,1)),lab0(natstore(k2,1),:),
     .            num0,sym,numb,nflag)
        if(nflag.eq.1)then
c        write(6,*)' k1,k2 ',k1,k2
c        write(6,*)' nat0(natstore(k2,1)) ',nat0(natstore(k2,1))
         do i=1,nat0(natstore(k2,1))
c         write(6,*)i,lab0(natstore(k2,1),i),num0(i)
         enddo
         stop
        endif

        netatoms=netatoms+nat0(natstore(k2,1))
        ne=ne-nch0(natstore(k2,1))
        do i=1,nat0(natstore(k2,1))
         ne=ne+num0(i)
        enddo

       else
        ncharge=ncharge+nch1(k2)
        call Atnum(60,nat1(k2),lab1(k2,:),
     .            num1,sym,numb,nflag)
        if(nflag.eq.1)then
c        write(6,*)' k1,k2 ',k1,k2
c        write(6,*)' nat1(k2) ',nat1(k2)
         do i=1,nat1(k2)
c         write(6,*)i,lab1(k2,i)
         enddo
         stop
        endif
        netatoms=netatoms+nat1(k2)
        ne=ne-nch1(k2)
        do i=1,nat1(k2)
         ne=ne+num1(i)
        enddo
       endif
c calculate the multiplicity
        ne2=ne/2
        ne3=ne-2*ne2
        if(ne3.eq.0)then
         multiplicity=1
        else
         multiplicity=2
        endif

c output the number of electrons for test information
       write(26,*)netatoms,ne,ca

       write(1,'(T1,A12,I7)')"!Isg_coeff= ",intsign
       write(1,*)

       write(1,*)ncharge,multiplicity

c write the coordinates
       if(numat(k1).eq.1)then
        do i=1,nat0(natstore(k1,1))
         write(1,100)lab0(natstore(k1,1),i),
     .               (c0(natstore(k1,1),i,k),k=1,3)
        enddo
       else
        do i=1,nat1(k1)
         write(1,100)lab1(k1,i),(c1(k1,i,k),k=1,3)
        enddo
       endif
c repeat for k2
       if(numat(k2).eq.1)then
        do i=1,nat0(natstore(k2,1))
         write(1,100)lab0(natstore(k2,1),i),
     .               (c0(natstore(k2,1),i,k),k=1,3)
        enddo
       else
        do i=1,nat1(k2)
         write(1,100)lab1(k2,i),(c1(k2,i,k),k=1,3)
        enddo
       endif

100    format(a2,3f13.6)
c add a blank line
       write(1,*)

c skip the embedded charges for this version of L1L1
       go to 3000

c add any point charges if necessary
      jc=0
      if(numcharges.gt.0)then
      do k=1,numcharges
c skip the charge group if it is contained in this L1 frag
        do j=1,numat(k1)
         if(natstore(k1,j).eq.numchggrps(k))go to 440
        enddo
        do j=1,numat(k2)
         if(natstore(k2,j).eq.numchggrps(k))go to 440
        enddo
c this charge group can be added 
       call filelabel(numchggrps(k),ca)
       n1=index(ca,' ')-1
       ca1='charge.'//ca(1:n1)//'.grp'
       ca2='charge.'//ca(1:n1)//'.grp'//'_atoms'
       open(unit=2,file=ca1,status='old')
       open(unit=7,file=ca2,status='old')
441    continue
       read(2,*,end=442)x1,y1,z1,ptch1
       write(1,101)x1,y1,z1,ptch1
       jc=jc+1
       read(7,*)num(jc),(nbg(jc,l),we(jc,l),l=1,num(jc))
       ptch(jc)=ptch1
       go to 441
442    continue
       close(unit=2)
       close(unit=7)
101    format(4f15.7)

c end the loop over charge files
440   continue
      enddo
c add a blank line
       write(1,*)
      endif
      close(unit=1)
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

      close(unit=4)
      return
      end
