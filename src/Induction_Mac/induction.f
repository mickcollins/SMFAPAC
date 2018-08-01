      program induction
      implicit double precision(a-h,o-z)
      parameter(natomf=120)

      integer, allocatable   :: nat0(:)
      integer, allocatable   :: numat(:), natstore(:,:)
      integer, allocatable   :: nfam(:), ifam(:,:)
      real*8, allocatable   :: d(:,:),cent(:,:),field(:,:)
      real*8, allocatable   :: c0(:,:,:),alpha(:,:,:)

      real*8, allocatable   :: fg(:,:,:,:)
      real*8, allocatable   :: g(:,:),h(:,:),w1(:,:,:)
      integer, allocatable   :: nat1(:),num1(:,:),natnum1(:,:,:)

      dimension ch(natomf),dip(natomf,3),coord(natomf,3)

      dimension t1(3),t2(3,3),t3(3,3,3),t4(3,3,3,3),t5(3,3,3,3,3)

      dimension q(3,3),o(3,3,3),hx(3,3,3,3),fd(3),fd1(3),dfield(3,3)

c the allowed maximum number of atoms in a group is natomf

c get the number of atoms
       open(unit=1,file='name.xyz',status='old')
       read(1,*)natom
       close(unit=1)
      open(unit=1,file='IN_JOBTYPE',status='old')
      read(1,*)
      read(1,*)jobtype
      close(unit=1)

      if(jobtype.gt.0)then
       allocate(g(natom,3))
       g=0
      endif
       if(jobtype.eq.2)then
        allocate(h(3*natom,3*natom))
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

c we need the number of original Level 1 frags
c and their signs and the groups they contain

      open(unit=1,file='OUT_FINAL_L1_DATA',status='old')
      read(1,*)
      read(1,*)
      read(1,*)nfrags
      read(1,*)
      allocate(numat(nfrags))
      allocate(natstore(nfrags,8))
      do n=1,nfrags
       read(1,*)i1,numat(n)
       read(1,*)(natstore(n,j),j=1,numat(n))
      enddo      

      close(unit=1)

c     write(6,*)nfrags
c      call readL1atoms(natom,NL1,nfrags1,nat1,num1,natnum1,w1,c,nsign)
c     write(6,*)nfrags1
c     stop
c get the Level0 group data

      open(unit=1,file='frags.out_Lev0',status='old')
      read(1,*)
      read(1,*)ngroups
      close(unit=1)

      allocate(c0(ngroups,60,3))
      allocate(alpha(ngroups,3,3))
      allocate(cent(ngroups,3))
      allocate(field(ngroups,3))
      allocate(nat0(ngroups))
      allocate(d(ngroups,3))
      allocate(nfam(ngroups))
      allocate(ifam(ngroups,60))

       call readgroupdata(ngroups,c0,nat0,alpha,cent)
       call readfamilies(ngroups,nfam,ifam)

c calculate the field at each centre
       do n=1,ngroups
        do k1=1,3
         field(n,k1)=0.d0
        enddo
       enddo

      if(jobtype.gt.0)then
       allocate(fg(ngroups,3,natom,3))
       fg=0.d0
      endif

       eps=0.005d0
c now read in the allowed interactions from the OUT_L1L1_data file

      open(unit=2,file='OUT_L1L1_data',status='unknown')
10    continue
      read(2,*,end=20)m1,m2,nsign


c m1 and m2 are the numbers of L1 frags

c     if(m1.ne.1.or.m2.ne.7)go to 10

c first use the moments of m2
      call readmoments(natomf,m2,ch,dip,coord,natom1)

c now calulate the field at each group in frag m1
      do l=1,natom1
      do ja=1,numat(m1)

      n=natstore(m1,ja)
      call Tfactor2(cent(n,:),coord(l,:),t0,t1,t2)
      call fld(1,t1,t2,t3,t4,t5,ch(l),dip(l,:),q,o,hx,fd)
      do k=1,3
       field(n,k)=field(n,k)+fd(k)*nsign
      enddo

      if(jobtype.gt.0)then

      do j=1,3
       cent(n,j)=cent(n,j)+eps
       call Tfactor2(cent(n,:),coord(l,:),t0,t1,t2)
       call fld(1,t1,t2,t3,t4,t5,ch(l),dip(l,:),q,o,hx,fd1)
       do k=1,3
       dfield(k,j)=(fd1(k)-fd(k))/eps
       enddo
       cent(n,j)=cent(n,j)-eps
      enddo

c allocate the gradients
       do k=1,3
       do j=1,3

       do i=1,nfam(n)
       fg(n,k,ifam(n,i),j)=fg(n,k,ifam(n,i),j)+nsign*dfield(k,j)/nfam(n)
       enddo

       do n1=1,num1(m2,l)
        fg(n,k,natnum1(m2,l,n1),j)=fg(n,k,natnum1(m2,l,n1),j)-
     .                         w1(m2,l,n1)*nsign*dfield(k,j)
       enddo

       enddo
       enddo
      endif

c end the l and ja loops
      enddo
      enddo


c repeat for the field at each group in frag m2
c use the moments of m1
      call readmoments(natomf,m1,ch,dip,coord,natom1)
      do l=1,natom1
      do ja=1,numat(m2)

      n=natstore(m2,ja)
      call Tfactor2(cent(n,:),coord(l,:),t0,t1,t2)
      call fld(1,t1,t2,t3,t4,t5,ch(l),dip(l,:),q,o,hx,fd)
      do k=1,3
       field(n,k)=field(n,k) + fd(k)*nsign
      enddo

      if(jobtype.gt.0)then

      do j=1,3
       cent(n,j)=cent(n,j)+eps
       call Tfactor2(cent(n,:),coord(l,:),t0,t1,t2)
       call fld(1,t1,t2,t3,t4,t5,ch(l),dip(l,:),q,o,hx,fd1)
       do k=1,3
       dfield(k,j)=(fd1(k)-fd(k))/eps
       enddo
       cent(n,j)=cent(n,j)-eps
      enddo

c allocate the gradients
       do k=1,3
       do j=1,3

       do i=1,nfam(n)
       fg(n,k,ifam(n,i),j)=fg(n,k,ifam(n,i),j)+nsign*dfield(k,j)/nfam(n)
       enddo
      
       do n1=1,num1(m1,l)
        fg(n,k,natnum1(m1,l,n1),j)=fg(n,k,natnum1(m1,l,n1),j)-
     .                         w1(m1,l,n1)*nsign*dfield(k,j)
       enddo

       enddo
       enddo
      endif

c end the l and ja loops
      enddo
      enddo


c finished with this L1..L1 term, so go to the next
      go to 10
20    continue
      close(unit=2)

c calc the perturbed dipole
       do n=1,ngroups
       do k1=1,3
        d(n,k1)=0.d0
        do k2=1,3
         d(n,k1)=d(n,k1)+alpha(n,k1,k2)*field(n,k2)
        enddo
       enddo
       enddo

c calc the induction energy
       Eind=0.d0
       do n=1,ngroups
       do k=1,3
       Eind=Eind-0.5*d(n,k)*field(n,k)
       enddo
       enddo

       open(unit=2,file='Induction_energy',status='unknown')
       write(2,*)Eind
       close(unit=2)

      open(unit=1,file='Induction_derivs',status='unknown')
      write(1,*)Eind
      write(1,*)jobtype

      if(jobtype.gt.0)then
c evaluate and write out the gradients
      do n=1,natom
      do k=1,3
      do m=1,ngroups
      do j1=1,3
      do j2=1,3
       g(n,k)=g(n,k)-0.5d0*(fg(m,j1,n,k)*alpha(m,j1,j2)*field(m,j2)
     .                    +field(m,j1)*alpha(m,j1,j2)*fg(m,j2,n,k))
      enddo
      enddo
      enddo
      enddo
      enddo

c write out
      if(jobtype.gt.0)then
      write(1,*)natom
      do n=1,natom
       write(1,101)(g(n,k),k=1,3)
      enddo
      endif
101   format(6f10.6)
      if(jobtype.eq.2)then
        n3=3*natom
        write(1,101)((h(i,j),j=1,i),i=1,n3)
      endif
      close(unit=1)

      endif
       end

