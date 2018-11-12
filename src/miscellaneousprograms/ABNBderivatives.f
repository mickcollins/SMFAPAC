      program ABNBderivatives

c this program evaluates the contribution to the energy derivatives
c wrt atomic coordinates that arises from the Level1-Level1
c ab initio jobs. The atoms are only those contained in the
c ab initio fragments.
c The input file(s) contain forces (hessians) from the Level 1
c fragment pair, followed by the corresponding forces (hessians)
c for each separate fragment.

      implicit double precision(a-h,o-z)

      real*8, allocatable   :: g(:,:),h(:,:),w1(:,:,:)
      integer, allocatable   :: nat1(:),num1(:,:),natnum1(:,:,:)

      integer, allocatable   :: num(:),natnum(:,:)
      real*8, allocatable   :: w(:,:)

      dimension abforce(3,200,3), fc(3,600,600)
c this assumes that two Level 1 fragments will not have a combined
c number of atoms exceeding 200. Very likely safe.

      character*20 ca1,ca2
      character*30 ca

c get the jobtype
      open(unit=1,file='IN_JOBTYPE',status='old')
      read(1,*)
      read(1,*)jobtype
      close(unit=1)

c this program should only be called if jobtype=1 (gradients)
c or jobtype=2 (hessians)
      if(jobtype.eq.0)then
       write(6,*)' called ABNBderivatives when jobtype  neither 1 or 2'
       stop
      endif

c get the number of atoms in the whole molecule
      open(unit=1,file='name.xyz',status='old')
      read(1,*)natom
      close(unit=1)

c allocate the gradient (g) and hessian (h) for the whole molecule
      allocate(g(natom,3))
      allocate(h(3*natom,3*natom))
      allocate(num(natom))
      allocate(natnum(natom,6))
      allocate(w(natom,6))

      do n=1,natom
      do k=1,3
       g(n,k)=0.d0
      enddo
      enddo

      h=0


c get the Lev1 identities
c OUT_Lev1_ATOMALLOCATION identifies the atom numbering in the whole
c molecule for the atoms in the Level 1 fragments
      open(unit=1,file='OUT_Lev1_ATOMALLOCATION',status='old')
      read(1,*,end=900)
      read(1,*)NL1
c NL1 is the number of Level 1 fragments
      read(1,*)
      read(1,*)
      allocate(nat1(NL1))
c nat1 = the number of atoms in the fragment
      allocate(num1(NL1,200))
c num1=1 for a real atom, num1=2 for a cap
c allocation assumes no more than 200 atoms in a Level 1 fragment
      allocate(natnum1(NL1,200,6))
c natnum1 = the whole molecule atom number
      allocate(w1(NL1,200,6))
c in the case of caps, w1 = the weight this cap has on a real atom

      do n=1,NL1
       read(1,*)i1,nat1(n)
       do m=1,nat1(n)
        read(1,*)num1(n,m),
     .           (natnum1(n,m,k),w1(n,m,k),k=1,num1(n,m))
       enddo
      enddo

      close(unit=1)

      open(unit=2,file='abforces',status='unknown')
c abforces contins the forces extracted from the nonbonded
c ab initio jobs
      if(jobtype.eq.2)then
      open(unit=3,file='abhessians',status='unknown')
c abhessians contains the hessions extracted from the nonbonded
c ab initio jobs

      endif 

c read in the number of ab jobs
c njobs is the total number of Level1-Level1 ab initio jobs
      read(2,*,end=127)njobs
      go to 128
127   continue
      stop
128   continue

       if(jobtype.eq.2)then
        read(3,*)
       endif

      do n=1,njobs
c read the two L1 indices
       read(2,*)n1,n2
c n1 and n2 identify the Level 1 fragments, whose atom numbers were
c read from OUT_Lev1_ATOMALLOCATION above

c read the sign
       read(2,*)nsign
c nsign is the coeffiient of the fragment pair in the sum
c over contributions to the nonbonded energy

c read the number of atoms in the ab file
       read(2,*)nat2
c nat2 = the number of atoms in the pair
c read the forces for the n1 n2 pair
c      read(2,*)
c      read(2,*)
c      read(2,*)
       do m=1,nat2
        read(2,*)i1,i2,(abforce(1,m,k),k=1,3)
       enddo
c i1 and i2 are not needed here
c      read(2,*)
c      read(2,*)
c      read(2,*)
c      read(2,*)
       read(2,*)
       read(2,*)
c read the forces for n1 
       do m=1,nat1(n1)
        read(2,*)i1,i2,(abforce(2,m,k),k=1,3)
       enddo
c i1 and i2 are not needed here
c      read(2,*)
c      read(2,*)
c      read(2,*)
c      read(2,*)
       read(2,*)
       read(2,*)
c read the forces for n2 
       do m=1,nat1(n2)
        read(2,*)i1,i2,(abforce(3,m,k),k=1,3)
       enddo
c i1 and i2 are not needed here
c      read(2,*)

       if(jobtype.eq.2)then
c read the second derivatives, lower triangle
c first the n1, n2 pair
        read(3,*)
        read(3,*)
        read(3,*)
        n3=3*nat2
       do i=1,n3
       do j=1,i
        read(3,*)fc(1,i,j)
       enddo
       enddo
c make fc symmetric
       do i=1,n3
       do j=i,n3
       fc(1,i,j)=fc(1,j,i)
       enddo
       enddo
c then read the n1 fragment data
        read(3,*)
        read(3,*)
        n3=3*nat1(n1)
       do i=1,n3
       do j=1,i
        read(3,*)fc(2,i,j)
       enddo
       enddo
c make fc symmetric
       do i=1,n3
       do j=i,n3
       fc(2,i,j)=fc(2,j,i)
       enddo
       enddo
        read(3,*)
        read(3,*)
c then read the n2 fragment data
        n3=3*nat1(n2)
       do i=1,n3
       do j=1,i
        read(3,*)fc(3,i,j)
       enddo
       enddo
c make fc symmetric
       do i=1,n3
       do j=i,n3
       fc(3,i,j)=fc(3,j,i)
       enddo
       enddo
c end the jojtype if
       endif

c allocate gradients
c first make a combined set of identities for n1 and n2
      num=0
      natnum=0
      w=0
      do m=1,nat1(n1)
       num(m)=num1(n1,m)
       do l=1,num(m)
       natnum(m,l)=natnum1(n1,m,l)
       w(m,l)=w1(n1,m,l)
       enddo
      enddo
      do m=1,nat1(n2)
       num(nat1(n1)+m)=num1(n2,m)
       do l=1,num1(n2,m)
        natnum(nat1(n1)+m,l)=natnum1(n2,m,l)
        w(nat1(n1)+m,l)=w1(n2,m,l)
       enddo
      enddo
      nat=nat1(n1)+nat1(n2)

       do m=1,nat1(n1)
       do j=1,num1(n1,m)
       do k=1,3
        g(natnum1(n1,m,j),k)=g(natnum1(n1,m,j),k)-nsign*
     .                   (abforce(1,m,k)-abforce(2,m,k))*w1(n1,m,j)
c the - sign in "-nsign" converts forces to gradients
c the weight w1=1 for real atoms, or a fraction for caps
       enddo
       enddo
       enddo
       do m=1,nat1(n2)
        m1=m+nat1(n1)
       do j=1,num1(n2,m)
       do k=1,3
        g(natnum1(n2,m,j),k)=g(natnum1(n2,m,j),k)-nsign*
     .                   (abforce(1,m1,k)-abforce(3,m,k))*w1(n2,m,j)
       enddo
       enddo
       enddo

c allocate hessians
c first from the pair
       if(jobtype.eq.2)then
        do k1=1,3
        do k2=1,3
        do i1=1,nat
        do i2=1,nat
         l1=3*(i1-1)+k1
         l2=3*(i2-1)+k2
        do j1=1,num(i1)
        do j2=1,num(i2)
         m1=3*(natnum(i1,j1)-1)+k1
         m2=3*(natnum(i2,j2)-1)+k2
         h(m1,m2)=h(m1,m2)+fc(1,l1,l2)*w(i1,j1)*w(i2,j2)*nsign
        enddo
        enddo
        enddo
        enddo
        enddo
        enddo
c then subtract the n1 fragment
        do k1=1,3
        do k2=1,3
        do i1=1,nat1(n1)
        do i2=1,nat1(n1)
         l1=3*(i1-1)+k1
         l2=3*(i2-1)+k2
        do j1=1,num1(n1,i1)
        do j2=1,num1(n1,i2)
         m1=3*(natnum1(n1,i1,j1)-1)+k1
         m2=3*(natnum1(n1,i2,j2)-1)+k2
         h(m1,m2)=h(m1,m2)-fc(2,l1,l2)*w1(n1,i1,j1)*w1(n1,i2,j2)*nsign
        enddo
        enddo
        enddo
        enddo
        enddo
        enddo
c then subtract the n2 fragment
        do k1=1,3
        do k2=1,3
        do i1=1,nat1(n2)
        do i2=1,nat1(n2)
         l1=3*(i1-1)+k1
         l2=3*(i2-1)+k2
        do j1=1,num1(n2,i1)
        do j2=1,num1(n2,i2)
         m1=3*(natnum1(n2,i1,j1)-1)+k1
         m2=3*(natnum1(n2,i2,j2)-1)+k2
         h(m1,m2)=h(m1,m2)-fc(3,l1,l2)*w1(n2,i1,j1)*w1(n2,i2,j2)*nsign
        enddo
        enddo
        enddo
        enddo
        enddo
        enddo

c end the jobtype if
      endif

c end the loop over ab jobs (njobs)
      enddo
      close(unit=2)

900   continue
c output
      open(unit=1,file='combinedABNBderivs',status='unknown')
      write(1,*)jobtype
      write(1,*)natom
      do n=1,natom
       write(1,100)(g(n,k),k=1,3)
      enddo
100   format(6f10.6)
      if(jobtype.eq.2)then
        n3=3*natom
        write(1,100)((h(i,j),j=1,i),i=1,n3)
      endif

      close(unit=1)
      end

