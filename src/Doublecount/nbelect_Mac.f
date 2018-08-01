      program L1L1elect
      implicit double precision(a-h,o-z)

c theis program operates just like Electro_Mac, except
c that it only calculates the electrostatic interaction
c between charged groups that are far apart.
c It reverses the sign of the interaction, so that these
c interactions are subtracted from the energy, as they are double
c counted in the FRAG part.



c reads in the cat OUT_L1L1_data file which contains lines like
c       k1     k2    nsign
c where k1 and k2 identify the L1 fragments that will be
c interacted electrostatically
c nsign is the sign coefficient

c We will read OUT_L1L1_data one line at a time
c we will open the
c nb.k1.0.cart
c nb.k2.0.cart
c files, and read in the multipole moments etc
c Then we simply call the elect subroutine to evalaute the
c multipole-multipole interactions

      real*8, allocatable    ::ch1(:),c1(:,:),d1(:,:),q1(:,:,:)
      real*8, allocatable    ::o1(:,:,:,:),h1(:,:,:,:,:)

      real*8, allocatable    ::ch2(:),c2(:,:),d2(:,:),q2(:,:,:)
      real*8, allocatable    ::o2(:,:,:,:),h2(:,:,:,:,:)

      real*8, allocatable    ::del1(:,:),del2(:,:)
      real*8, allocatable    ::secd1(:,:,:,:),secd2(:,:,:,:),
     .                         secd3(:,:,:,:)

      character*20 ca, ca1

      integer, allocatable  ::nat0(:)
      integer, allocatable  ::numchggrps(:)
      real*8, allocatable   ::c0(:,:,:)
      integer, allocatable  ::numat(:),natstore(:,:)
      character*2, lab0

      real*8, allocatable   :: g(:,:),h(:,:),w1(:,:,:)
      integer, allocatable   :: nat1(:),num1(:,:),natnum1(:,:,:)

      dimension totk(3)

c     integer, allocatable   :: num(:),natnum(:,:)
c     real*8, allocatable   :: w(:,:)

c get the jobtype
      open(unit=1,file='IN_JOBTYPE',status='old')
      read(1,*)
      read(1,*)jobtype
      close(unit=1)

c the IN_COMBINE file specifies whether electrostatic forces are
c used
c     open(unit=1,file='IN_COMBINE',status='old')
c     read(1,*)
c     read(1,*)nforce
c     close(unit=1)
c DO NOT USE GRADIENTS
c     if(nforce.lt.3)jobtype=0


      tote=0

      numcharges=0
c we will need to know which groups contain charges
      open(unit=1,file='OUT_CHARGEDGROUPS',status='unknown')
      read(1,*,end=333)
      read(1,*)numcharges
      if(numcharges.gt.0)then
       allocate(numchggrps(numcharges))
       read(1,*)
       do k=1,numcharges
        read(1,*)numchggrps(k)
        read(1,*)
        read(1,*)
        read(1,*)
       enddo
      endif
333   continue
      close(unit=1)

c if there are less than 2 charges, there is nothing to do
      if(numcharges.lt.2)go to 999


c we need the number of groups
      open(unit=1,file='frags.out_Lev0',status='old')
      read(1,*)
      read(1,*)ngroups
      allocate(nat0(ngroups))
      allocate(c0(ngroups,60,3))
      close(unit=1)
c if there are charges, we will need the coordinates of the groups
c which contain charges


c read in the Level 0 (group) coordinates
      do n=1,ngroups
       call filelabel(n,ca)
       ca1='Lev0_COORD'//ca
       open(unit=1,file=ca1,status='old')
       read(1,*)
       m=1
15     continue
       read(1,*,end=25)lab0,(c0(n,m,k),k=1,3)
       m=m+1
       go to 15
25     continue
       nat0(n)=m-1
       close(unit=1)
      enddo
100   format(a2,3f13.6)
c we will need the L1 fragment group numbers
      open(unit=8,file='OUT_FINAL_L1_DATA',status='unknown')
      read(8,*)
      read(8,*)
      read(8,*)nL1frags
      allocate(numat(nL1frags))
      allocate(natstore(nL1frags,8))

      read(8,*)
      do n=1,nL1frags
       read(8,*)n1,numat(n)
       read(8,*)(natstore(n,j),j=1,numat(n))
      enddo
      close(unit=8)

c if jobtype > 0, we need data to allocate gradients and
c second derivatives
      if(jobtype.gt.0)then
c get the number of atoms
       open(unit=1,file='name.xyz',status='old')
       read(1,*)natom
       close(unit=1)

       allocate(g(natom,3)) 
       g=0
      endif
       if(jobtype.eq.2)then
        allocate(h(3*natom,3*natom))
c       allocate(num(natom))
c       allocate(natnum(natom,6))
c       allocate(w(natom,6))
        h=0
       endif

c get the Lev1 identities
       open(unit=1,file='OUT_Lev1_ATOMALLOCATION',status='old')
       read(1,*)
       read(1,*)NL1
       read(1,*)
       read(1,*)

       allocate(nat1(NL1))
       allocate(num1(NL1,200))
       allocate(natnum1(NL1,200,6))
       allocate(w1(NL1,200,6))

       do n=1,NL1
        read(1,*)i1,nat1(n)
        do m=1,nat1(n)
         read(1,*)num1(n,m),
     .           (natnum1(n,m,k),w1(n,m,k),k=1,num1(n,m))
        enddo
       enddo

       close(unit=1)

       do n=1,NL1
       do m=1,nat1(n)
        if(num1(n,m).eq.1)w1(n,m,1)=1.d0
        if(num1(n,m).eq.2)w1(n,m,2)=1.d0-w1(n,m,1)
       enddo
       enddo
c now we read in the electrostatic interactions L1..L1, one at a time

      open(unit=2,file='OUT_CHARGE_CHARGE',status='unknown')
c OUT_CHARGE_CHARGE replaces OUT_L1L1_data in Electro_Mac

c begin the line by line loop
1     continue
      read(2,*,end=2)m1,m2,nsign
c     write(6,*)m1,m2,nsign

c make the file name, first k1
       call filelabel(m1,ca1)
       n1=index(ca1,' ')-1
       ca='nb.'//ca1(1:n1)//'.0'//'.cart'
       open(unit=1,file=ca,status='old')
       read(1,*)
c read the number of atoms
       read(1,*)natom1
        if(natom1.ne.nat1(m1))then
         write(6,*)' problem with natom1 ',m1,nat1(m1),natom1
         stop
        endif
       allocate(ch1(natom1))
       allocate(c1(natom1,3))
       allocate(d1(natom1,3))
       allocate(q1(natom1,3,3))
       allocate(o1(natom1,3,3,3))
       allocate(h1(natom1,3,3,3,3))
       allocate(del1(natom1,3))

c start reading the data
      do n=1,natom1
       read(1,*)
       do k=1,3
        read(1,*)c1(n,k)
       enddo
       read(1,*)
       read(1,*)ch1(n)
       read(1,*)
       do k=1,3
        read(1,*)d1(n,k)
       enddo
       read(1,*)
       do k1=1,3
       do k2=k1,3
       read(1,*)q1(n,k1,k2)
       enddo
       enddo
       read(1,*)
       do k1=1,3
       do k2=k1,3
       do k3=k2,3
       read(1,*)o1(n,k1,k2,k3)
       enddo
       enddo
       enddo
       read(1,*)
       do k1=1,3
       do k2=k1,3
       do k3=k2,3
       do k4=k3,3
       read(1,*)h1(n,k1,k2,k3,k4)
       enddo
       enddo
       enddo
       enddo
       call symm(q1(n,:,:),o1(n,:,:,:),h1(n,:,:,:,:))
c end the n loop
      enddo

      close(unit=1)
c In Electro_Mac
c  if there are charges, we need to see if this L1 frag contains
c  group with charges. The electrostatic interaction of charged
c  groups will all other groups has already been done, so the
c  moments of such groups must be set to zero
c here
c  the charge-neutral interactions should also not be done
c  but we zero the neutral groups

      if(numcharges.gt.0)then
       
       do j=1,numat(m1)
        nmatch=0
        do k=1,numcharges
         if(numchggrps(k).eq.natstore(m1,j))nmatch=1
        enddo
        if(nmatch.eq.0)then
c       write(6,*)m1,natstore(m1,j)
         call zero(natom1,c1,ch1,d1,q1,o1,h1,
     .             nat0(natstore(m1,j)),c0(natstore(m1,j),:,:))
        endif
       enddo
      endif
c make the file name, now k2
       call filelabel(m2,ca1)
       n1=index(ca1,' ')-1
       ca='nb.'//ca1(1:n1)//'.0'//'.cart'
       open(unit=1,file=ca,status='old')

       read(1,*)
c read the number of atoms
       read(1,*)natom2
        if(natom2.ne.nat1(m2))then
         write(6,*)' problem with natom2 ',m2,nat1(m2),natom2
         stop
        endif
       allocate(ch2(natom2))
       allocate(c2(natom2,3))
       allocate(d2(natom2,3))
       allocate(q2(natom2,3,3))
       allocate(o2(natom2,3,3,3))
       allocate(h2(natom2,3,3,3,3))
       allocate(del2(natom2,3))


c start reading the data
      do n=1,natom2
       read(1,*)
       do k=1,3
        read(1,*)c2(n,k)
       enddo
       read(1,*)
       read(1,*)ch2(n)
       read(1,*)
       do k=1,3
        read(1,*)d2(n,k)
       enddo
       read(1,*)
       do k1=1,3
       do k2=k1,3
       read(1,*)q2(n,k1,k2)
       enddo
       enddo
       read(1,*)
       do k1=1,3
       do k2=k1,3
       do k3=k2,3
       read(1,*)o2(n,k1,k2,k3)
       enddo
       enddo
       enddo
       read(1,*)
       do k1=1,3
       do k2=k1,3
       do k3=k2,3
       do k4=k3,3
       read(1,*)h2(n,k1,k2,k3,k4)
       enddo
       enddo
       enddo
       enddo

       call symm(q2(n,:,:),o2(n,:,:,:),h2(n,:,:,:,:))
c end the n loop
      enddo

      close(unit=1)
c if there are charges, we need to see if this L1 frag contains
c group with charges. The electrostatic interaction of charged
c groups will all other groups has already been done, so the
c moments of such groups must be set to zero

c here the neutral groups are zeroed
      if(numcharges.gt.0)then

       do j=1,numat(m2)
        nmatch=0
       do k=1,numcharges
        if(numchggrps(k).eq.natstore(m2,j))nmatch=1
       enddo
        if(nmatch.eq.0)then
         call zero(natom2,c2,ch2,d2,q2,o2,h2,
     .             nat0(natstore(m2,j)),c0(natstore(m2,j),:,:))
        endif      
       enddo
      endif
c       write(6,*)m2,natstore(m2,j)
c now we have all the multipoles

c convert coordinates to Bohr
      Bohr=1.88972598860
      do n=1,natom1
      do k=1,3
       c1(n,k)=c1(n,k)*Bohr
      enddo
      enddo
      do n=1,natom2
      do k=1,3
       c2(n,k)=c2(n,k)*Bohr
      enddo
      enddo

c we use only up to quadrupoles for the fast code
c the multipole tensors do not have to be symmetrized

c simply call subroutine elect
c jobtype=0, energies only
c jobtype=1, also gradients
c jobtype=2 second derivatives


        allocate(secd1(natom1,3,natom1,3))
        allocate(secd2(natom1,3,natom2,3))
        allocate(secd3(natom2,3,natom2,3))

      call elect(natom1,c1,ch1,d1,q1,o1,h1,
     .           natom2,c2,ch2,d2,q2,o2,h2,
     .           elen,del1,del2,secd1,secd2,secd3,jobtype)
      write(66,*)m1,m2,nsign,elen

c     do n=1,natom1
c      write(50,*)(del1(n,k),k=1,3)
c     enddo
c      write(50,*)
c     do n=1,natom2
c      write(50,*)(del2(n,k),k=1,3)
c     enddo

      tote=tote+elen*nsign

      if(jobtype.gt.0)then

c make the gradients sum to zero
      do k=1,3
       totk(k)=0.d0
      enddo

c allocate the gradients
       do k=1,3
       do m=1,nat1(m1)
        totk(k)=totk(k)+del1(m,k)
       do j=1,num1(m1,m)
        g(natnum1(m1,m,j),k)=g(natnum1(m1,m,j),k)+
     .                      del1(m,k)*w1(m1,m,j)*nsign
c       totk(k)=totk(k)+del1(m,k)*w1(m1,m,j)
       enddo
       enddo
       enddo

       do k=1,3
       do m=1,nat1(m2)
        totk(k)=totk(k)+del2(m,k)
       do j=1,num1(m2,m)
        g(natnum1(m2,m,j),k)=g(natnum1(m2,m,j),k)+
     .                      del2(m,k)*w1(m2,m,j)*nsign
c       totk(k)=totk(k)+del2(m,k)*w1(m2,m,j)
       enddo
       enddo
       enddo

c      write(34,*)(totk(k),k=1,3)
c end the jobtype if
      endif 

      if(jobtype.eq.2)then

        do k1=1,3
        do k2=1,3
        do i1=1,natom1
c       do i2=1,natom1
        i2=i1
        do j1=1,num1(m1,i1)
        do j2=1,num1(m1,i2)
         l1=3*(natnum1(m1,i1,j1)-1)+k1
         l2=3*(natnum1(m1,i2,j2)-1)+k2
         h(l1,l2)=h(l1,l2)+secd1(i1,k1,i2,k2)*w1(m1,i1,j1)*
     .                     w1(m1,i2,j2)*nsign
        enddo
        enddo
c       enddo
        enddo
        enddo
        enddo

        do k1=1,3
        do k2=1,3
        do i1=1,natom2
c       do i2=1,natom2
        i2=i1
        do j1=1,num1(m2,i1)
        do j2=1,num1(m2,i2)
         l1=3*(natnum1(m2,i1,j1)-1)+k1
         l2=3*(natnum1(m2,i2,j2)-1)+k2
         h(l1,l2)=h(l1,l2)+secd3(i1,k1,i2,k2)*w1(m2,i1,j1)*
     .                     w1(m2,i2,j2)*nsign
        enddo
        enddo
c       enddo
        enddo
        enddo
        enddo

        do k1=1,3
        do k2=1,3
        do i1=1,natom1
        do i2=1,natom2
        do j1=1,num1(m1,i1)
        do j2=1,num1(m2,i2)
         l1=3*(natnum1(m1,i1,j1)-1)+k1
         l2=3*(natnum1(m2,i2,j2)-1)+k2
         h(l1,l2)=h(l1,l2)+secd2(i1,k1,i2,k2)*w1(m1,i1,j1)*
     .                     w1(m2,i2,j2)*nsign
         h(l2,l1)=h(l2,l1)+secd2(i1,k1,i2,k2)*w1(m1,i1,j1)*
     .                     w1(m2,i2,j2)*nsign
        enddo
        enddo
        enddo
        enddo
        enddo
        enddo

c end the jobtype=2 if
      endif

c deallocate the arrays
       deallocate(ch1)
       deallocate(c1)
       deallocate(d1)
       deallocate(q1)
       deallocate(o1)
       deallocate(h1)
       deallocate(del1)

       deallocate(ch2)
       deallocate(c2)
       deallocate(d2)
       deallocate(q2)
       deallocate(o2)
       deallocate(h2)
       deallocate(del2)

       deallocate(secd1)
       deallocate(secd2)
       deallocate(secd3)

      go to 1
2     continue
      close(unit=2)

c check sums
c      do n1=1,natom
c      do k1=1,3
c      do k2=1,3
c       sumk=0.d0
c       do n2=1,natom
c        i1=3*(n1-1)+k1
c        i2=3*(n2-1)+k2
c        sumk=sumk+h(i1,i2)
c       enddo
c       write(34,*)sumk
c      enddo
c      enddo
c      enddo

c output
      if(jobtype.gt.0)then

c remove centre of mass forces
      do k=1,3
       sumk=0.d0
       do n=1,natom
        sumk=sumk+g(n,k)
       enddo
c      write(6,*)sumk,natom
       sumk=sumk/dble(float(natom))
       do n=1,natom
        g(n,k)=g(n,k)-sumk
       enddo
      enddo

      endif

      open(unit=1,file='Doublecount_derivs',status='unknown')
      write(1,*)-tote
c sign reversed to remove double counting
      write(1,*)jobtype
      if(jobtype.gt.0)then
      write(1,*)natom
      do n=1,natom
       write(1,*)(-g(n,k),k=1,3)
c sign reversed to remove double counting
      enddo
      endif
101   format(6f10.6)
      if(jobtype.eq.2)then
        n3=3*natom
        write(1,101)((-h(i,j),j=1,i),i=1,n3)
c sign reversed to remove double counting
c      do i=1,n3
c      do j=1,i
c       write(1,*)h(i,j)
c      enddo
c      enddo
      endif
      close(unit=1)

999   continue
      end
      




