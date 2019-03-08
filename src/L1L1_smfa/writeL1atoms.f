      subroutine writeL1atoms(nffinal1,ngroups,nfextra,numat,natstore)

      implicit double precision(a-h,o-z)

      dimension numat(2*nffinal1),natstore(2*nffinal1,8)
      dimension numat1(2*nffinal1),natstore1(2*nffinal1,8)

      dimension natL1(nffinal1),natatomL1(nffinal1,60),ncapL1(nffinal1),
     .          icapL1(nffinal1,24,2)
      dimension fcapL1(nffinal1,24)

      dimension natL0(ngroups),natatomL0(ngroups,60),ncapL0(ngroups),
     .          icapL0(ngroups,24,2)
      dimension fcapL0(ngroups,24)





c  the purpose of this subroutine is to create an output file that contains
c the data needed to allocate the force on the kth atom in the nth L1
c fragment to one or more atoms in the original molecule. For most L1
c fragment atoms, there is a 1-to-1 correspondence, but for caps
c it is 1-to-2.

c get the original L1 data, but we will have to cope with
c cancelations done in subroutine cancelL1

      open(unit=8,file='OUT_FRAGMENTS_Lev1',status='old')
      read(8,*)nffinal
      do i=1,nffinal
      read(8,*)numat1(i)
c     read(8,*)(natstore1(i,k),k=1,numat(i))
c corrected this line on 27/09/18
      read(8,*)(natstore1(i,k),k=1,numat1(i))
      enddo
      close(unit=8)

c sort for later comparison with natstore
      do i=1,nffinal
       call piksrt(numat1(i),natstore1(i,:))
      enddo

      open(unit=3,file='frags.out_Lev1',status='old')
      read(3,*)
      read(3,*)nf
      if(nf.ne.nffinal)then
       write(6,*)' OUT_FRAGMENTS_Lev1 and frags.out_Lev1 do not match'
       stop
      endif
      read(3,*)
      read(3,*)
      read(3,*)
      do m=1,nffinal
      read(3,*)natL1(m)
      enddo
      read(3,*)
      do m=1,nffinal
      read(3,*)(natatomL1(m,i),i=1,natL1(m))
      enddo
c read the cap data
      do m=1,nffinal
       ncapL1(m)=0
      enddo
      read(3,*)
20    read(3,*,end=21)m1,n1,n2,f
      ncapL1(m1)=ncapL1(m1)+1
      icapL1(m1,ncapL1(m1),1)=n1
      icapL1(m1,ncapL1(m1),2)=n2
      fcapL1(m1,ncapL1(m1))=f
      go to 20
21    continue
      close(unit=3)

c get the original L0 data

      open(unit=3,file='frags.out_Lev0',status='old')
      read(3,*)
      read(3,*)ngroups0
      read(3,*)
      read(3,*)
      read(3,*)
      do m=1,ngroups
      read(3,*)natL0(m)
      enddo
      read(3,*)
      do m=1,ngroups
      read(3,*)(natatomL0(m,i),i=1,natL0(m))
      enddo
c read the cap data
      do m=1,ngroups
       ncapL0(m)=0
      enddo
      read(3,*)
22    read(3,*,end=23)m1,n1,n2,f
      ncapL0(m1)=ncapL0(m1)+1
      icapL0(m1,ncapL0(m1),1)=n1
      icapL0(m1,ncapL0(m1),2)=n2
      fcapL0(m1,ncapL0(m1))=f
      go to 22
23    continue
      close(unit=3)

c each real atom in a final L1 fragment has one corresponding real atom
c caps have 2

c write out all the necessary data
#ifdef __GFORTRAN__
      open(unit=1,file='OUT_Lev1_ATOMALLOCATION',status='unknown')
#else
      open(unit=1,file='OUT_Lev1_ATOMALLOCATION',status='unknown',
     .            buffered='YES')
#endif
      write(1,*)' The number of final L1 frags from L1L1_mac.f'
      write(1,*)nfextra
      write(1,*)' For each of ',nfextra,' fragments:the number of atoms'
      write(1,*)' followed by the allocation of each atom'
      do n=1,nfextra
       match=0
       if(numat(n).gt.1)then
c identify the original Lev1 frag
        do j=1,nffinal
         if(numat1(j).ne.numat(n))go to 24
         do k=1,numat(n)
          if(natstore1(j,k).ne.natstore(n,k))go to 24
         enddo
c we have a match for the original Lev1 frag j
         match=1
c write out the atoms for original fragment j
         write(1,*)n,natL1(j)+ncapL1(j)
         do k=1,natL1(j)
         write(1,100)1,natatomL1(j,k),1.d0
         enddo
         do k=1,ncapL1(j)
         write(1,101)2,icapL1(j,k,1),1.d0-fcapL1(j,k),icapL1(j,k,2),fcapL1(j,k)
         enddo
         go to 25
c end the nffinal loop
24     enddo
25     continue
       if(match.eq.0)then
        write(6,*)' no match for multi group L1 frag'
        stop
       endif
c end the numat if
       endif
       if(numat(n).eq.1)then
c the original Lev0 group is used
        j=natstore(n,1)
c write out the atoms for original fragment j
         write(1,*)n,natL0(j)+ncapL0(j)
         do k=1,natL0(j)
         write(1,100)1,natatomL0(j,k),1.d0
         enddo
         do k=1,ncapL0(j)
         write(1,101)2,icapL0(j,k,1),1.d0-fcapL0(j,k),icapL0(j,k,2),fcapL0(j,k)
         enddo
c end the second numat if
       endif
c end the nfextra loop
      enddo
100    format(i4,i8,f12.6)
101    format(i4,i8,f12.6,i8,f12.6)

      return
      end
