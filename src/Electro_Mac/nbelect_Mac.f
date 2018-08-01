      program L1L1elect
      implicit double precision(a-h,o-z)

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

      real*8, allocatable    ::ch1(:,:),c1(:,:,:),d1(:,:,:),q1(:,:,:,:)
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

      Bohr=1.88972598860

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

c find max atoms in a level 1 fragment
      L1max=0
       do n=1,NL1
        if(nat1(n).gt.L1max)L1max=nat1(n)
       enddo

       allocate(ch1(NL1,L1max))
       allocate(c1(NL1,L1max,3))
       allocate(d1(NL1,L1max,3))
       allocate(q1(NL1,L1max,3,3))
c read in the multipoles
       do m=1,NL1
       call filelabel(m,ca1)
       n1=index(ca1,' ')-1
       ca='nb.'//ca1(1:n1)//'.0'//'.cart'
       open(unit=1,file=ca,status='old')
       read(1,*)
       read(1,*)natom1
        if(natom1.ne.nat1(m))then
         write(6,*)' problem with natom1 ',m,nat1(m),natom1
         stop
        endif
      do n=1,natom1
       read(1,*)
       do k=1,3
        read(1,*)c1(m,n,k)
       enddo
       read(1,*)
       read(1,*)ch1(m,n)
       read(1,*)
       do k=1,3
        read(1,*)d1(m,n,k)
       enddo
       read(1,*)
       do k1=1,3
       do k2=k1,3
       read(1,*)q1(m,n,k1,k2)
       q1(m,n,k2,k1)=q1(m,n,k1,k2)
       enddo
       enddo
       read(1,*)
       do k1=1,3
       do k2=k1,3
       do k3=k2,3
        read(1,*)
       enddo
       enddo
       enddo
       read(1,*)
       do k1=1,3
       do k2=k1,3
       do k3=k2,3
       do k4=k3,3
        read(1,*)
       enddo
       enddo
       enddo
       enddo

      enddo
       
      close(unit=1)
c if there are charges, we need to see if this L1 frag contains
c group with charges. The electrostatic interaction of charged
c groups will all other groups has already been done, so the
c moments of such groups must be set to zero
      if(numcharges.gt.0)then
       do k=1,numcharges
       do j=1,numat(m)
        if(numchggrps(k).eq.natstore(m,j))then
       call zero1(L1max,natom1,c1(m,:,:),ch1(m,:),d1(m,:,:),q1(m,:,:,:),
     .             nat0(natstore(m,j)),c0(natstore(m,j),:,:))
        endif
       enddo
       enddo
      endif
c convert coordinates to Bohr
      do n=1,nat1(m)
      do k=1,3
       c1(m,n,k)=c1(m,n,k)*Bohr
      enddo
      enddo

c close the m loop
      enddo

c now we read in the electrostatic interactions L1..L1, one at a time

      open(unit=2,file='OUT_L1L1_data',status='unknown')

c begin the line by line loop
1     continue
      read(2,*,end=2)m1,m2,nsign

c we use only up to quadrupoles for the fast code
c the multipole tensors do not have to be symmetrized

c simply call subroutine elect
c jobtype=0, energies only
c jobtype=1, also gradients
c jobtype=2 second derivatives

        allocate(del1(nat1(m1),3))
        allocate(del2(nat1(m2),3))
        allocate(secd1(nat1(m1),3,nat1(m1),3))
        allocate(secd2(nat1(m1),3,nat1(m2),3))
        allocate(secd3(nat1(m2),3,nat1(m2),3))

      call elect(L1max,nat1(m1),c1(m1,:,:),ch1(m1,:),d1(m1,:,:),
     .           q1(m1,:,:,:),
     .           nat1(m2),c1(m2,:,:),ch1(m2,:),d1(m2,:,:),q1(m2,:,:,:),
     .           elen,del1,del2,secd1,secd2,secd3,jobtype)

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
        do i1=1,nat1(m1)
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
        do i1=1,nat1(m2)
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
        do i1=1,nat1(m1)
        do i2=1,nat1(m2)
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
       deallocate(del1)
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

      open(unit=1,file='Electrostatic_derivs',status='unknown')
      write(1,*)tote
      write(1,*)jobtype
      if(jobtype.gt.0)then
      write(1,*)natom
      do n=1,natom
       write(1,*)(g(n,k),k=1,3)
      enddo
      endif
101   format(6f10.6)
      if(jobtype.eq.2)then
        n3=3*natom
        write(1,101)((h(i,j),j=1,i),i=1,n3)
c      do i=1,n3
c      do j=1,i
c       write(1,*)h(i,j)
c      enddo
c      enddo
      endif
      close(unit=1)


      end
      




