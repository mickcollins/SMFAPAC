
      program Disp
c This program calculates the dispersion energy for those
c long range group=group interactions that were not
c included as ab initio interactions in either the
c main fragments or the ab initio Level1-Level1 terms

      implicit double precision(a-h,o-z)

      integer, allocatable   :: natbig(:),numbig(:,:)
      integer, allocatable   ::nom(:,:)

      real*8, allocatable   ::c0(:,:,:),paa(:),alpha(:,:,:)
      real*8, allocatable   :: g(:,:),h(:,:),w0(:,:,:)
      real*8, allocatable   :: del1(:,:),del2(:,:)
      real*8, allocatable   :: ssecd1(:,:),ssecd2(:,:)
      real*8, allocatable   :: ssecd3(:,:)

      integer, allocatable   :: n0(:),nat0(:),num0(:,:),natnum0(:,:,:)


      character*2, allocatable  ::lab0(:,:)

      character*20 ca
      character*16 ca1

      dimension natstore(16)

      bigtol=10.d0

c get the number of atoms in the whole molecule
      open(unit=1,file='name.xyz',status='old')
      read(1,*)natom
      close(unit=1)

c allocate the gradient and hessian arrays
      allocate(g(natom,3))
      g=0.d0
      allocate(h(3*natom,3*natom))
      h=0.d0

c get the jobtype
      open(unit=1,file='IN_JOBTYPE',status='old')
      read(1,*)
      read(1,*)jobtype
      close(unit=1)
c jobtype=0 for energy only, jobtype=1 for grdient,
c jobtype=2 for hessian

c get the Level0 group data
      open(unit=1,file='frags.out_Lev0',status='old')
c Level 0 fragmentation produced individual groups
      read(1,*)
      read(1,*)ngroups
c ngroup = the number of groups in the molecule

      natomf=60
c natomf is the maximum allowed number of atoms in a group

      allocate(n0(ngroups))
      allocate(nat0(ngroups))
      allocate(c0(ngroups,natomf,3))
      allocate(lab0(ngroups,natomf))
      allocate(paa(ngroups))
      allocate(alpha(ngroups,3,3))
      allocate(nom(ngroups,ngroups))
      allocate(num0(ngroups,natomf))
      allocate(w0(ngroups,natomf,6))
      allocate(natnum0(ngroups,natomf,6))

      allocate(del1(natomf,3))
      allocate(del2(natomf,3))
      allocate(ssecd1(3*natomf,3*natomf))
      allocate(ssecd2(3*natomf,3*natomf))
      allocate(ssecd3(3*natomf,3*natomf))


      read(1,*)
      read(1,*)
      read(1,*)
      do n=1,ngroups
       read(1,*)n0(n)
c n0 is the number of real atoms in the group
      enddo
      read(1,*)
      do n=1,ngroups
       do m=1,n0(n)
        num0(n,m)=1
        w0(n,m,1)=1.d0
       enddo
       read(1,*)(natnum0(n,m,1),m=1,n0(n))
c natnum0 = the atom number in the whole molecule
      enddo
      read(1,*)

c read in the caps
50    continue
      read(1,*,end=51)ng1,m1,m2,f1
c ng1= the group number
c m1 and m1 are the atoms which define the cap, m1 is in
c group ng1
c f1 is the weight
      n0(ng1)=n0(ng1)+1
c the number of atoms in group ng1 increases by 1 for this cap
      num0(ng1,n0(ng1))=2
      natnum0(ng1,n0(ng1),1)=m1
      w0(ng1,n0(ng1),1)=1.d0-f1
      natnum0(ng1,n0(ng1),2)=m2
      w0(ng1,n0(ng1),2)=f1
      go to 50
51    continue

      close(unit=1)

c read in the Level 0 (group) coordinates
      do n=1,ngroups
       call filelabel(n,ca)
       ca1='Lev0_COORD'//ca
       open(unit=1,file=ca1,status='old')
       read(1,*)nch0
       m=1
1      continue
c      read(1,100,end=2)lab0(n,m),(c0(n,m,k),k=1,3)
       read(1,*,end=2)lab0(n,m),(c0(n,m,k),k=1,3)
       m=m+1
       go to 1
2      continue
       nat0(n)=m-1
c nat0 should be the number of real atoms + caps
       close(unit=1)
      enddo
100   format(a2,3f13.6)

c convert to Bohr
      bohr=1.d0/1.8897259886d0
      do n=1,ngroups
      do m=1,nat0(n)
      do k=1,3
       c0(n,m,k)=c0(n,m,k)/bohr
      enddo
      enddo
      enddo

c check n0 and nat0 are correct
      do n=1,ngroups
       if(n0(n).ne.nat0(n))then
        write(6,*)' confusion about the number of atoms in group ',n
        stop
       endif
      enddo

c read in the group polarizabilities and PAAA factors
      open(unit=1,file='Polarisabilities.out',status='old')
      do n=1,ngroups
       read(1,*)ii,alpha(n,1,1),alpha(n,2,1),alpha(n,2,2),
     .             alpha(n,3,1),alpha(n,3,2),alpha(n,3,3)
       alpha(n,1,2)=alpha(n,2,1)
       alpha(n,1,3)=alpha(n,3,1)
       alpha(n,2,3)=alpha(n,3,2)
      enddo
      close(unit=1)

      open(unit=1,file='PAA.out',status='old')
      do n=1,ngroups
       read(1,*)ii,paa(n)
      enddo
      close(unit=1)

c see Addicoat and Collins for how the zero frequency polarizability
c and P factor for each group are used.

c the dispersion interaction is calculated as an interaction between
c groups. The interaction between groups that have both been involved
c in ab initio calculations have already had their interaction
c accounted for. So, we find these pairs of groups and "exclude" them

c read in the Level 3 fragment group numbers
      open(unit=9,file='OUT_FRAGMENTS_LevX',status='old')
      read(9,*)nfbig
c OUT_FRAGMENTS_Lev3 contains data for the main fragments
c nfbig = the number of main fragments
      allocate(natbig(nfbig))
      allocate(numbig(nfbig,natomf))

      do i=1,nfbig
      read(9,*)natbig(i)
c natbig = the number of groups in the fragment
      read(9,*)(numbig(i,k),k=1,natbig(i))
c numbig = the group numbers
      enddo
      close(unit=9)


c form exclusion matrix
c nom(n1,n2)=1 means the dispersion interaction between groups
c n1 and n2 will be calculated
      do n1=1,ngroups
      do n2=1,ngroups
      nom(n1,n2)=1
      enddo
      enddo

      do i=1,nfbig
       do n1=1,natbig(i)
       do n2=1,natbig(i)
        nom(numbig(i,n1),numbig(i,n2))=0
       enddo
       enddo
      enddo
c get the ab initio nonbonded group numbers
c these will be added to the exclusion matrix
      open(unit=4,file='OUT_ABGROUPS',status='unknown')
c OUT_ABGROUPS lists the groups involved in nonbonded ab initio
c calculations
11    continue
      read(4,*,end=12)numabgroups
c numabgroups= the number of groups in a single nonbonded ab intio calc
      read(4,*)(natstore(k),k=1,numabgroups)
c natstore lists the group numbers

c adjust the exclusion matrix
      do j1=1,numabgroups-1
      do j2=j1+1,numabgroups
       nom(natstore(j1),natstore(j2))=0
       nom(natstore(j2),natstore(j1))=0
      enddo
      enddo
      go to 11
12    continue
      close(unit=4)


c now do the dispersion interactions
c loop over pairs of groups, using the exclusion atrix to skip some
c pairs
      energy=0.d0
      do n1=1,ngroups-1
      do n2=n1+1,ngroups

      if(nom(n1,n2).eq.0)go to 2000

c subroutine dispe calculates the dispersion interaction
c and its derivatives for apair of groups


      call dispe(natomf,nat0(n1),c0(n1,:,:),nat0(n2),c0(n2,:,:),paa(n1),
     .          paa(n2),
     .          alpha(n1,:,:),alpha(n2,:,:),ed,del1,del2,
     .          ssecd1,ssecd2,ssecd3,jobtype)


      do i=1,3*nat0(n1)
      do j=1,3*nat0(n1)
       if(abs(ssecd1(i,j)).gt.1.d2)then
        write(6,*)ssecd1(i,j)
        stop
       endif
      enddo
      enddo
      do i=1,3*nat0(n1)
      do j=1,3*nat0(n2)
       if(abs(ssecd2(i,j)).gt.1.d2)then
        write(6,*)ssecd1(i,j)
        stop
       endif
      enddo
      enddo
      do i=1,3*nat0(n2)
      do j=1,3*nat0(n2)
       if(abs(ssecd1(i,j)).gt.1.d2)then
        write(6,*)ssecd3(i,j)
        stop
       endif
      enddo
      enddo


c allocate gradients and seconds
      energy=energy+ed
c gradients
c for real atoms, num0=1 and w0=1
c for caps, num0=2 and the weight w0 is a fraction

      if(jobtype.ge.1)then
      do m=1,nat0(n1)
      do j=1,num0(n1,m)
      do k=1,3
       g(natnum0(n1,m,j),k)=del1(m,k)*w0(n1,m,j)
      enddo
      enddo
      enddo
      do m=1,nat0(n2)
      do j=1,num0(n2,m)
      do k=1,3
       g(natnum0(n2,m,j),k)=del1(m,k)*w0(n2,m,j)
      enddo
      enddo
      enddo
      endif
c second derivatives
      if(jobtype.eq.2)then

        do k1=1,3
        do k2=1,3
        do i1=1,nat0(n1)
        do i2=1,nat0(n1)
         m1=3*(i1-1)+k1
         m2=3*(i2-1)+k2
        do j1=1,num0(n1,i1)
        do j2=1,num0(n1,i2)
         l1=3*(natnum0(n1,i1,j1)-1)+k1
         l2=3*(natnum0(n1,i2,j2)-1)+k2
         h(l1,l2)=h(l1,l2)+ssecd1(m1,m2)*w0(n1,i1,j1)*
     .                     w0(n1,i2,j2)
        enddo
        enddo
        enddo
        enddo
        enddo
        enddo

        do k1=1,3
        do k2=1,3
        do i1=1,nat0(n2)
        do i2=1,nat0(n2)
         m1=3*(i1-1)+k1
         m2=3*(i2-1)+k2
        do j1=1,num0(n2,i1)
        do j2=1,num0(n2,i2)
         l1=3*(natnum0(n2,i1,j1)-1)+k1
         l2=3*(natnum0(n2,i2,j2)-1)+k2
         h(l1,l2)=h(l1,l2)+ssecd3(m1,m2)*w0(n2,i1,j1)*
     .                     w0(n2,i2,j2)
        enddo
        enddo
        enddo
        enddo
        enddo
        enddo

        do k1=1,3
        do k2=1,3
        do i1=1,nat0(n1)
        do i2=1,nat0(n2)
         m1=3*(i1-1)+k1
         m2=3*(i2-1)+k2
        do j1=1,num0(n1,i1)
        do j2=1,num0(n2,i2)
         l1=3*(natnum0(n1,i1,j1)-1)+k1
         l2=3*(natnum0(n2,i2,j2)-1)+k2
         h(l1,l2)=h(l1,l2)+ssecd2(m1,m2)*w0(n1,i1,j1)*
     .                     w0(n2,i2,j2)
         h(l2,l1)=h(l2,l1)+ssecd2(m1,m2)*w0(n1,i1,j1)*
     .                     w0(n2,i2,j2)
        enddo
        enddo
        enddo
        enddo
        enddo
        enddo


c end the jobtype if
      endif
2000  continue


c end the loops over n1,n2
      enddo
      enddo

c  output
      open(unit=1,file='Dispderivatives',status='unknown')
      write(1,*)energy
      write(1,*)jobtype
      write(1,*)natom
      if(jobtype.ge.1)then
      do m=1,natom
       write(1,101)(g(m,k),k=1,3)
      enddo
      endif
      if(jobtype.eq.2)then
        n3=3*natom
        write(1,101)((h(i,j),j=1,i),i=1,n3)
101   format(6f10.6)
      endif
      close(unit=1)

      end


      subroutine dispe(natomf,natom1,x1,natom2,x2,paa1,paa2,
     .                alpha1,alpha2,elen,del1,del2,
     .                ssecd1,ssecd2,ssecd3,jobtype)
      implicit double precision(a-h,o-z)

      dimension del1(natomf,3),del2(natomf,3)
      dimension x1(natomf,3),x2(natomf,3),y1(3),y2(3)

      dimension den(3),d2en(3,3)
      dimension ssecT(6*natomf,6*natomf)
      dimension ssecd1(3*natomf,3*natomf),ssecd2(3*natomf,3*natomf),
     .          ssecd3(3*natomf,3*natomf)   


      dimension alpha1(3,3),alpha2(3,3)

      dimension t1(3),t2(3,3),t3(3,3,3),t4(3,3,3,3)

      dimension x(2*natomf,3),Pro(6*natomf,6*natomf),gra(6*natomf)
     .       ,grb(6*natomf),tmp(6*natomf,6*natomf)

c calculate the centroid of the two groups and store in y1 and y2

      do k=1,3
       y1(k)=0.d0
      do n=1,natom1
       y1(k)=y1(k)+x1(n,k)/natom1
      enddo
      enddo
      do k=1,3
       y2(k)=0.d0
      do n=1,natom2
       y2(k)=y2(k)+x2(n,k)/natom2
      enddo
      enddo
      
      call Tfactor4(y1,y2,t0,t1,t2,t3,t4)
c Tfactor4 calculates the inverse distance of y1 nd y2 and
c derivatives of this (t0) up to 4th order
 
      f=-2.d0*paa1*paa2/(paa1+paa2)
c see Addicoat and Collins for this ratio of P factors

      elen=0.d0
      do k1=1,3
       den(k1)=0.d0
      do k2=1,3
       d2en(k1,k2)=0.d0
      enddo
      enddo
c elen, den and d2en are the dispersion energy and its derivatives

      ssecd1=0
      ssecd2=0
      ssecd3=0

      do k1=1,3
      do k2=1,3
      do k3=1,3
      do k4=1,3
       ans1=alpha1(k1,k3)*alpha2(k2,k4)*f
       e1=t2(k1,k2)*t2(k3,k4)*ans1
       elen=elen+e1
c see Addicoat and Collins for this product of polarisabilities
c and derivatives of t0
      do k5=1,3
       de1=(t3(k5,k1,k2)*t2(k3,k4)+t2(k1,k2)*t3(k5,k3,k4))
       den(k5)=den(k5)+de1*ans1
      do k6=1,3
       e2=t4(k6,k5,k1,k2)*t2(k3,k4)+t3(k5,k1,k2)*t3(k6,k3,k4)+
     .    t3(k6,k1,k2)*t3(k5,k3,k4)+t2(k1,k2)*t4(k6,k5,k3,k4)
       d2en(k5,k6)=d2en(k5,k6)+e2*ans1
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo

c     write(10,*)1.d0/t0,elen


      if(jobtype.gt.0)then   

      do k=1,3
      do n=1,natom1
       del1(n,k)=den(k)/natom1
      enddo
      do n=1,natom2
       del2(n,k)=-den(k)/natom2
      enddo
      enddo

      endif

      if(jobtype.eq.2)then

c all second derivatives
      do k1=1,3
      do k2=1,3

       do j1=1,natom1
        i1=3*(j1-1)+k1
       do j2=1,natom1
        i2=3*(j2-1)+k2
        ssecd1(i1,i2)=d2en(k1,k2)/natom1**2
       enddo
       enddo

       do j1=1,natom1
        i1=3*(j1-1)+k1
       do j2=1,natom2
        i2=3*(j2-1)+k2
        ssecd2(i1,i2)=-d2en(k1,k2)/(natom1*natom2)
       enddo
       enddo

       do j1=1,natom2
        i1=3*(j1-1)+k1
       do j2=1,natom2
        i2=3*(j2-1)+k2
        ssecd3(i1,i2)=d2en(k1,k2)/natom2**2
       enddo
       enddo

      enddo
      enddo

      endif


      return
      end

      subroutine Tfactor4(x1,x2,t0,t1,t2,t3,t4)
      implicit double precision(a-h,o-z)

      dimension x1(3),x2(3),r(3),dx(3)
      dimension t1(3),t2(3,3),t3(3,3,3),t4(3,3,3,3)
c     .          ,t5(3,3,3,3,3),t6(3,3,3,3,3,3)

c x1 and x2 are the coordinates of the two atoms

c r is the vector from x2 to x1
c and r0 is the distance, |r|.

      r0=0.d0
      do k=1,3
      r(k)=x1(k)-x2(k)
      r0=r0+r(k)**2
      enddo
      r0=sqrt(r0)
 
c  First we equate the t factors, t1, t2, etc to the part ofeach
c that has no delta functions

      t0=1.d0/r0

      do k1=1,3
      t1(k1)=-t0*r(k1)/r0**2
      do k2=1,3
      t2(k1,k2)=-3.d0*t1(k1)*r(k2)/r0**2
      do k3=1,3
      t3(k1,k2,k3)=-5.d0*t2(k1,k2)*r(k3)/r0**2
      do k4=1,3
      t4(k1,k2,k3,k4)=-7.d0*t3(k1,k2,k3)*r(k4)/r0**2
c     do k5=1,3
c     t5(k1,k2,k3,k4,k5)=-9.d0*t4(k1,k2,k3,k4)*r(k5)/r0**2
c     do k6=1,3
c     t6(k1,k2,k3,k4,k5,k6)=-11.d0*t5(k1,k2,k3,k4,k5)*r(k6)/r0**2
c     enddo
c     enddo
      enddo
      enddo
      enddo
      enddo

c Now we add on the parts that have one delta function, by using the
c same index for two array indices.

      do k1=1,3
      t2(k1,k1)=t2(k1,k1)-1.d0/r0**3
      do k2=1,3
      t3(k1,k1,k2)=t3(k1,k1,k2)+3.d0*r(k2)/r0**5
      t3(k1,k2,k1)=t3(k1,k2,k1)+3.d0*r(k2)/r0**5
      t3(k2,k1,k1)=t3(k2,k1,k1)+3.d0*r(k2)/r0**5
      do k3=1,3
      t4(k1,k1,k2,k3)=t4(k1,k1,k2,k3)-15.d0*r(k2)*r(k3)/r0**7
      t4(k1,k2,k1,k3)=t4(k1,k2,k1,k3)-15.d0*r(k2)*r(k3)/r0**7
      t4(k1,k2,k3,k1)=t4(k1,k2,k3,k1)-15.d0*r(k2)*r(k3)/r0**7
      t4(k2,k1,k1,k3)=t4(k2,k1,k1,k3)-15.d0*r(k2)*r(k3)/r0**7
      t4(k2,k1,k3,k1)=t4(k2,k1,k3,k1)-15.d0*r(k2)*r(k3)/r0**7
      t4(k2,k3,k1,k1)=t4(k2,k3,k1,k1)-15.d0*r(k2)*r(k3)/r0**7
c     do k4=1,3
c     t5(k1,k1,k2,k3,k4)=
c    .t5(k1,k1,k2,k3,k4)+105.d0*r(k2)*r(k3)*r(k4)/r0**9
c     t5(k1,k2,k1,k3,k4)=
c    .t5(k1,k2,k1,k3,k4)+105.d0*r(k2)*r(k3)*r(k4)/r0**9
c     t5(k1,k2,k3,k1,k4)=
c    .t5(k1,k2,k3,k1,k4)+105.d0*r(k2)*r(k3)*r(k4)/r0**9
c     t5(k1,k2,k3,k4,k1)=
c    .t5(k1,k2,k3,k4,k1)+105.d0*r(k2)*r(k3)*r(k4)/r0**9
c     t5(k2,k1,k1,k3,k4)=
c    .t5(k2,k1,k1,k3,k4)+105.d0*r(k2)*r(k3)*r(k4)/r0**9
c     t5(k2,k1,k3,k1,k4)=
c    .t5(k2,k1,k3,k1,k4)+105.d0*r(k2)*r(k3)*r(k4)/r0**9
c     t5(k2,k1,k3,k4,k1)=
c    .t5(k2,k1,k3,k4,k1)+105.d0*r(k2)*r(k3)*r(k4)/r0**9
c     t5(k2,k3,k1,k1,k4)=
c    .t5(k2,k3,k1,k1,k4)+105.d0*r(k2)*r(k3)*r(k4)/r0**9
c     t5(k2,k3,k1,k4,k1)=
c    .t5(k2,k3,k1,k4,k1)+105.d0*r(k2)*r(k3)*r(k4)/r0**9
c     t5(k2,k3,k4,k1,k1)=
c    .t5(k2,k3,k4,k1,k1)+105.d0*r(k2)*r(k3)*r(k4)/r0**9
c     do k5=1,3
c     t6(k1,k1,k2,k3,k4,k5)=
c    .t6(k1,k1,k2,k3,k4,k5)-945.d0*r(k2)*r(k3)*r(k4)*r(k5)/r0**11
c     t6(k1,k2,k1,k3,k4,k5)=
c    .t6(k1,k2,k1,k3,k4,k5)-945.d0*r(k2)*r(k3)*r(k4)*r(k5)/r0**11
c     t6(k1,k2,k3,k1,k4,k5)=
c    .t6(k1,k2,k3,k1,k4,k5)-945.d0*r(k2)*r(k3)*r(k4)*r(k5)/r0**11
c     t6(k1,k2,k3,k4,k1,k5)=
c    .t6(k1,k2,k3,k4,k1,k5)-945.d0*r(k2)*r(k3)*r(k4)*r(k5)/r0**11
c     t6(k1,k2,k3,k4,k5,k1)=
c    .t6(k1,k2,k3,k4,k5,k1)-945.d0*r(k2)*r(k3)*r(k4)*r(k5)/r0**11
c     t6(k2,k1,k1,k3,k4,k5)=
c    .t6(k2,k1,k1,k3,k4,k5)-945.d0*r(k2)*r(k3)*r(k4)*r(k5)/r0**11
c     t6(k2,k1,k3,k1,k4,k5)=
c    .t6(k2,k1,k3,k1,k4,k5)-945.d0*r(k2)*r(k3)*r(k4)*r(k5)/r0**11
c     t6(k2,k1,k3,k4,k1,k5)=
c    .t6(k2,k1,k3,k4,k1,k5)-945.d0*r(k2)*r(k3)*r(k4)*r(k5)/r0**11
c     t6(k2,k1,k3,k4,k5,k1)=
c    .t6(k2,k1,k3,k4,k5,k1)-945.d0*r(k2)*r(k3)*r(k4)*r(k5)/r0**11
c     t6(k2,k3,k1,k1,k4,k5)=
c    .t6(k2,k3,k1,k1,k4,k5)-945.d0*r(k2)*r(k3)*r(k4)*r(k5)/r0**11
c     t6(k2,k3,k1,k4,k1,k5)=
c    .t6(k2,k3,k1,k4,k1,k5)-945.d0*r(k2)*r(k3)*r(k4)*r(k5)/r0**11
c     t6(k2,k3,k1,k4,k5,k1)=
c    .t6(k2,k3,k1,k4,k5,k1)-945.d0*r(k2)*r(k3)*r(k4)*r(k5)/r0**11
c     t6(k2,k3,k4,k1,k1,k5)=
c    .t6(k2,k3,k4,k1,k1,k5)-945.d0*r(k2)*r(k3)*r(k4)*r(k5)/r0**11
c     t6(k2,k3,k4,k1,k5,k1)=
c    .t6(k2,k3,k4,k1,k5,k1)-945.d0*r(k2)*r(k3)*r(k4)*r(k5)/r0**11
c     t6(k2,k3,k4,k5,k1,k1)=
c    .t6(k2,k3,k4,k5,k1,k1)-945.d0*r(k2)*r(k3)*r(k4)*r(k5)/r0**11
c     enddo
c     enddo
      enddo
      enddo
      enddo

c  Now we add on the parts that have two delta functions in a similar
c  way.

      do k1=1,3
      do k2=1,3
      t4(k1,k1,k2,k2)=t4(k1,k1,k2,k2)+3.d0/r0**5
      t4(k1,k2,k1,k2)=t4(k1,k2,k1,k2)+3.d0/r0**5
      t4(k1,k2,k2,k1)=t4(k1,k2,k2,k1)+3.d0/r0**5
c     do k3=1,3
c     t5(k1,k1,k2,k2,k3)=
c    .t5(k1,k1,k2,k2,k3)-15.d0*r(k3)/r0**7
c     t5(k1,k1,k2,k3,k2)=
c    .t5(k1,k1,k2,k3,k2)-15.d0*r(k3)/r0**7
c     t5(k1,k1,k3,k2,k2)=
c    .t5(k1,k1,k3,k2,k2)-15.d0*r(k3)/r0**7

c     t5(k1,k2,k1,k2,k3)=
c    .t5(k1,k2,k1,k2,k3)-15.d0*r(k3)/r0**7
c     t5(k1,k2,k1,k3,k2)=
c    .t5(k1,k2,k1,k3,k2)-15.d0*r(k3)/r0**7
c     t5(k1,k3,k1,k2,k2)=
c    .t5(k1,k3,k1,k2,k2)-15.d0*r(k3)/r0**7

c     t5(k1,k2,k2,k1,k3)=
c    .t5(k1,k2,k2,k1,k3)-15.d0*r(k3)/r0**7
c     t5(k1,k2,k3,k1,k2)=
c    .t5(k1,k2,k3,k1,k2)-15.d0*r(k3)/r0**7
c     t5(k1,k3,k2,k1,k2)=
c    .t5(k1,k3,k2,k1,k2)-15.d0*r(k3)/r0**7

c     t5(k1,k2,k2,k3,k1)=
c    .t5(k1,k2,k2,k3,k1)-15.d0*r(k3)/r0**7
c     t5(k1,k2,k3,k2,k1)=
c    .t5(k1,k2,k3,k2,k1)-15.d0*r(k3)/r0**7
c     t5(k1,k3,k2,k2,k1)=
c    .t5(k1,k3,k2,k2,k1)-15.d0*r(k3)/r0**7

c     t5(k3,k1,k1,k2,k2)=
c    .t5(k3,k1,k1,k2,k2)-15.d0*r(k3)/r0**7
c     t5(k3,k1,k2,k1,k2)=
c    .t5(k3,k1,k2,k1,k2)-15.d0*r(k3)/r0**7
c     t5(k3,k1,k2,k2,k1)=
c    .t5(k3,k1,k2,k2,k1)-15.d0*r(k3)/r0**7
c     do k4=1,3
c     t6(k1,k1,k2,k2,k3,k4)=
c    .t6(k1,k1,k2,k2,k3,k4)+105.d0*r(k3)*r(k4)/r0**9
c     t6(k1,k1,k2,k3,k2,k4)=
c    .t6(k1,k1,k2,k3,k2,k4)+105.d0*r(k3)*r(k4)/r0**9
c     t6(k1,k1,k2,k3,k4,k2)=
c    .t6(k1,k1,k2,k3,k4,k2)+105.d0*r(k3)*r(k4)/r0**9
c     t6(k1,k1,k3,k2,k2,k4)=
c    .t6(k1,k1,k3,k2,k2,k4)+105.d0*r(k3)*r(k4)/r0**9
c     t6(k1,k1,k3,k2,k4,k2)=
c    .t6(k1,k1,k3,k2,k4,k2)+105.d0*r(k3)*r(k4)/r0**9
c     t6(k1,k1,k3,k4,k2,k2)=
c    .t6(k1,k1,k3,k4,k2,k2)+105.d0*r(k3)*r(k4)/r0**9

c     t6(k1,k2,k1,k2,k3,k4)=
c    .t6(k1,k2,k1,k2,k3,k4)+105.d0*r(k3)*r(k4)/r0**9
c     t6(k1,k2,k1,k3,k2,k4)=
c    .t6(k1,k2,k1,k3,k2,k4)+105.d0*r(k3)*r(k4)/r0**9
c     t6(k1,k2,k1,k3,k4,k2)=
c    .t6(k1,k2,k1,k3,k4,k2)+105.d0*r(k3)*r(k4)/r0**9
c     t6(k1,k3,k1,k2,k2,k4)=
c    .t6(k1,k3,k1,k2,k2,k4)+105.d0*r(k3)*r(k4)/r0**9
c     t6(k1,k3,k1,k2,k4,k2)=
c    .t6(k1,k3,k1,k2,k4,k2)+105.d0*r(k3)*r(k4)/r0**9
c     t6(k1,k3,k1,k4,k2,k2)=
c    .t6(k1,k3,k1,k4,k2,k2)+105.d0*r(k3)*r(k4)/r0**9

c     t6(k1,k2,k2,k1,k3,k4)=
c    .t6(k1,k2,k2,k1,k3,k4)+105.d0*r(k3)*r(k4)/r0**9
c     t6(k1,k2,k3,k1,k2,k4)=
c    .t6(k1,k2,k3,k1,k2,k4)+105.d0*r(k3)*r(k4)/r0**9
c     t6(k1,k2,k3,k1,k4,k2)=
c    .t6(k1,k2,k3,k1,k4,k2)+105.d0*r(k3)*r(k4)/r0**9
c     t6(k1,k3,k2,k1,k2,k4)=
c    .t6(k1,k3,k2,k1,k2,k4)+105.d0*r(k3)*r(k4)/r0**9
c     t6(k1,k3,k2,k1,k4,k2)=
c    .t6(k1,k3,k2,k1,k4,k2)+105.d0*r(k3)*r(k4)/r0**9
c     t6(k1,k3,k4,k1,k2,k2)=
c    .t6(k1,k3,k4,k1,k2,k2)+105.d0*r(k3)*r(k4)/r0**9

c     t6(k1,k2,k2,k3,k1,k4)=
c    .t6(k1,k2,k2,k3,k1,k4)+105.d0*r(k3)*r(k4)/r0**9
c     t6(k1,k2,k3,k2,k1,k4)=
c    .t6(k1,k2,k3,k2,k1,k4)+105.d0*r(k3)*r(k4)/r0**9
c     t6(k1,k2,k3,k4,k1,k2)=
c    .t6(k1,k2,k3,k4,k1,k2)+105.d0*r(k3)*r(k4)/r0**9
c     t6(k1,k3,k2,k2,k1,k4)=
c    .t6(k1,k3,k2,k2,k1,k4)+105.d0*r(k3)*r(k4)/r0**9
c     t6(k1,k3,k2,k4,k1,k2)=
c    .t6(k1,k3,k2,k4,k1,k2)+105.d0*r(k3)*r(k4)/r0**9
c     t6(k1,k3,k4,k2,k1,k2)=
c    .t6(k1,k3,k4,k2,k1,k2)+105.d0*r(k3)*r(k4)/r0**9

c     t6(k1,k2,k2,k3,k4,k1)=
c    .t6(k1,k2,k2,k3,k4,k1)+105.d0*r(k3)*r(k4)/r0**9
c     t6(k1,k2,k3,k2,k4,k1)=
c    .t6(k1,k2,k3,k2,k4,k1)+105.d0*r(k3)*r(k4)/r0**9
c     t6(k1,k2,k3,k4,k2,k1)=
c    .t6(k1,k2,k3,k4,k2,k1)+105.d0*r(k3)*r(k4)/r0**9
c     t6(k1,k3,k2,k2,k4,k1)=
c    .t6(k1,k3,k2,k2,k4,k1)+105.d0*r(k3)*r(k4)/r0**9
c     t6(k1,k3,k2,k4,k2,k1)=
c    .t6(k1,k3,k2,k4,k2,k1)+105.d0*r(k3)*r(k4)/r0**9
c     t6(k1,k3,k4,k2,k2,k1)=
c    .t6(k1,k3,k4,k2,k2,k1)+105.d0*r(k3)*r(k4)/r0**9
c new code 6/105

c     t6(k3,k1,k1,k2,k2,k4)=
c    .t6(k3,k1,k1,k2,k2,k4)+105.d0*r(k3)*r(k4)/r0**9
c     t6(k3,k1,k1,k2,k4,k2)=
c    .t6(k3,k1,k1,k2,k4,k2)+105.d0*r(k3)*r(k4)/r0**9
c     t6(k3,k1,k1,k4,k2,k2)=
c    .t6(k3,k1,k1,k4,k2,k2)+105.d0*r(k3)*r(k4)/r0**9

c     t6(k3,k1,k2,k1,k2,k4)=
c    .t6(k3,k1,k2,k1,k2,k4)+105.d0*r(k3)*r(k4)/r0**9
c     t6(k3,k1,k2,k1,k4,k2)=
c    .t6(k3,k1,k2,k1,k4,k2)+105.d0*r(k3)*r(k4)/r0**9
c     t6(k3,k1,k4,k1,k2,k2)=
c    .t6(k3,k1,k4,k1,k2,k2)+105.d0*r(k3)*r(k4)/r0**9

c     t6(k3,k1,k2,k2,k1,k4)=
c    .t6(k3,k1,k2,k2,k1,k4)+105.d0*r(k3)*r(k4)/r0**9
c     t6(k3,k1,k2,k4,k1,k2)=
c    .t6(k3,k1,k2,k4,k1,k2)+105.d0*r(k3)*r(k4)/r0**9
c     t6(k3,k1,k4,k2,k1,k2)=
c    .t6(k3,k1,k4,k2,k1,k2)+105.d0*r(k3)*r(k4)/r0**9

c     t6(k3,k1,k2,k2,k4,k1)=
c    .t6(k3,k1,k2,k2,k4,k1)+105.d0*r(k3)*r(k4)/r0**9
c     t6(k3,k1,k2,k4,k2,k1)=
c    .t6(k3,k1,k2,k4,k2,k1)+105.d0*r(k3)*r(k4)/r0**9
c     t6(k3,k1,k4,k2,k2,k1)=
c    .t6(k3,k1,k4,k2,k2,k1)+105.d0*r(k3)*r(k4)/r0**9

c     t6(k3,k4,k1,k1,k2,k2)=
c    .t6(k3,k4,k1,k1,k2,k2)+105.d0*r(k3)*r(k4)/r0**9
c     t6(k3,k4,k1,k2,k1,k2)=
c    .t6(k3,k4,k1,k2,k1,k2)+105.d0*r(k3)*r(k4)/r0**9
c     t6(k3,k4,k1,k2,k2,k1)=
c    .t6(k3,k4,k1,k2,k2,k1)+105.d0*r(k3)*r(k4)/r0**9

c     enddo
c     enddo
      enddo
      enddo

c  t6 has a part with three delta functions

c     do k1=1,3
c     do k2=1,3
c     do k3=1,3
c     t6(k1,k1,k2,k2,k3,k3)=
c    .t6(k1,k1,k2,k2,k3,k3)-15.d0/r0**7
c     t6(k1,k2,k1,k2,k3,k3)=
c    .t6(k1,k2,k1,k2,k3,k3)-15.d0/r0**7
c     t6(k1,k2,k2,k1,k3,k3)=
c    .t6(k1,k2,k2,k1,k3,k3)-15.d0/r0**7
c     t6(k1,k2,k2,k3,k1,k3)=
c    .t6(k1,k2,k2,k3,k1,k3)-15.d0/r0**7
c     t6(k1,k2,k2,k3,k3,k1)=
c    .t6(k1,k2,k2,k3,k3,k1)-15.d0/r0**7

c     t6(k1,k1,k2,k3,k2,k3)=
c    .t6(k1,k1,k2,k3,k2,k3)-15.d0/r0**7
c     t6(k1,k1,k2,k3,k3,k2)=
c    .t6(k1,k1,k2,k3,k3,k2)-15.d0/r0**7

c     t6(k1,k2,k1,k3,k2,k3)=
c    .t6(k1,k2,k1,k3,k2,k3)-15.d0/r0**7
c     t6(k1,k2,k1,k3,k3,k2)=
c    .t6(k1,k2,k1,k3,k3,k2)-15.d0/r0**7

c     t6(k1,k2,k3,k1,k2,k3)=
c    .t6(k1,k2,k3,k1,k2,k3)-15.d0/r0**7
c     t6(k1,k2,k3,k1,k3,k2)=
c    .t6(k1,k2,k3,k1,k3,k2)-15.d0/r0**7

c     t6(k1,k2,k3,k2,k1,k3)=
c    .t6(k1,k2,k3,k2,k1,k3)-15.d0/r0**7
c     t6(k1,k2,k3,k3,k1,k2)=
c    .t6(k1,k2,k3,k3,k1,k2)-15.d0/r0**7


c     t6(k1,k2,k3,k2,k3,k1)=
c    .t6(k1,k2,k3,k2,k3,k1)-15.d0/r0**7
c     t6(k1,k2,k3,k3,k2,k1)=
c    .t6(k1,k2,k3,k3,k2,k1)-15.d0/r0**7

c     enddo
c     enddo
c     enddo

      return
      end



      subroutine filelabel(k,ca)
      implicit double precision(a-h,o-z)

c this subroutines returns a character variable ca
c for the integer k

      character*20 ca
      character*1 I0(0:10),anum(10),ca1,ca2,ca3,ca4

      I0(10)='0'
      I0(0)='0'
      I0(1)='1'
      I0(2)='2'
      I0(3)='3'
      I0(4)='4'
      I0(5)='5'
      I0(6)='6'
      I0(7)='7'
      I0(8)='8'
      I0(9)='9'

      if(k.le.9)then
       ca=I0(k)
      endif
      if(k.gt.9.and.k.le.99)then
      k1=k/10
      k2=k-k1*10
      if(k2.eq.0)k2=10
      ca1=I0(k1)
      ca2=I0(k2)
      ca=ca1//ca2
      endif
      if(k.gt.99.and.k.le.999)then
      k1=k/100
      k3=k/10 - INT(k1)*10
      k4=k - INT(k1)*100 - INT(k3)*10
      ca1=I0(k1)
      ca2=I0(k3)
      ca3=I0(k4)
      ca=ca1//ca2//ca3
      endif
      if(k.gt.999.and.k.le.9999)then
      k1=k/1000
      k3=k/100 - INT(k1)*10
      k4=k/10 - INT(k1)*100 - INT(k3)*10
      k5=k-INT(k1)*1000-k3*100-k4*10
      ca1=I0(k1)
      ca2=I0(k3)
      ca3=I0(k4)
      ca4=I0(k5)
      ca=ca1//ca2//ca3//ca4

      endif

      return
      end

