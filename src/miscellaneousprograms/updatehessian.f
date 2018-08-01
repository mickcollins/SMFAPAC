      program updatehessian
      implicit double precision(a-h,o-z)

      real*8, allocatable   :: g(:),h(:,:),dx(:),delh(:,:)
      real*8, allocatable   :: gold(:),hold(:,:),z(:),delg(:)
      real*8, allocatable   :: hsr1(:,:),hbfgs(:,:),hpsb(:,:),h2(:,:)
      real*8, allocatable   :: amas(:),coord(:,:)

      character*2, allocatable   :: lab(:)
      character*20 step, filename

      character*50 com1,com2,com3,com4,com5,com6,com7,com8

      open(unit=1,file='xyzFILENAME',status='old')
      read(1,*)
      read(1,*)natom
      close(unit=1)

      n3=3*natom

      allocate(z(n3))
      allocate(g(n3))
      allocate(h(n3,n3))
      allocate(gold(n3))
      allocate(hold(n3,n3))
      allocate(hsr1(n3,n3))
      allocate(hbfgs(n3,n3))
      allocate(hpsb(n3,n3))
      allocate(h2(n3,n3))
      allocate(delg(n3))
      allocate(delh(n3,n3))
      allocate(dx(n3))
      allocate(coord(3,natom))
      allocate(lab(natom))
      allocate(amas(natom))

c The optimisation script saves files containing
c the energy, gradient and heesian at each step in
c the optimisation. For step "0", the hessian is calculated
c but for nstep > 0, it is an updated hessian.

c read in whether to search for a minimum or a TS
      open(unit=2,file='IN_JOBTYPE_original',status='old')
      read(2,*)
      read(2,*)Nsearch
      close(unit=2)
c Nsearch=3 for opt, 4 for TS

c read in the number of negative eigenvalues, as this determines
c which update formula to use
      open(unit=8,file='NEGEIGVALS',status='old')
      read(8,*)negeigs
      close(unit=8)

c The file STEPNUMBER contains the current step number in the
c optimisation 

      open(unit=1,file='STEPNUMBER',status='old')
      read(1,*)nstep
      close(unit=1)
c set nstep to the previous step
      nstep=nstep-1
c get the previous energy, gradient and hessian
      call filelabel(nstep,step)
      n1=index(step,' ')-1

      filename="energy."//step(1:n1)
      open(unit=1,file=filename,status='old')
      read(1,*)energy
      close(unit=1)

      filename="gradient."//step(1:n1)
      open(unit=1,file=filename,status='old')
      do n=1,natom
      do k=1,3
       nk=3*(n-1)+k
       read(1,*)gold(nk)
      enddo
      enddo
      close(unit=1)

      filename="hessian."//step(1:n1)
      open(unit=1,file=filename,status='old')
      read(1,100)((hold(i,j),j=1,i),i=1,n3)
      close(unit=1)
100   format(6f10.6)
      filename="deltax."//step(1:n1)
      open(unit=1,file=filename,status='old')
      do n=1,natom
      do k=1,3
       nk=3*(n-1)+k
       read(1,*)dx(nk)
      enddo
      enddo

c sym sec derivs
      do i=1,n3
      do j=1,i
       hold(j,i)=hold(i,j)
      enddo
      enddo

c read combinedderivs
      open(unit=1,file='combinedderivs',status='old')
      read(1,50)com1
      read(1,*)jobtype
      read(1,50)com2
      read(1,*)
      read(1,50)com3
      do n=1,natom
      read(1,*)amas(n)
      enddo
      read(1,50)com4
      do n=1,natom
      read(1,51)lab(n)
      enddo
50    format(a50)
51    format(a2)
      read(1,50)com5
      read(1,*)energynew
      read(1,50)com6
      do n=1,natom
      read(1,33)(coord(k,n),k=1,3)
      enddo
33      format(1x,3D25.16)
      read(1,50)com7
      do n=1,natom
      do k=1,3
       nk=3*(n-1)+k
       read(1,33)g(nk)
      enddo
      enddo
      read(1,50)com8
      close(unit=1)

c the combinedderivs file has a zero hessian, so we don't read it in

c if nstep=1 originally, now nstep=0
c we do not update the hessian as the old gradients are HF
c and the new gradients may be something else
      if(nstep.eq.0)then
       do i=1,n3
       do j=1,n3
        h(i,j)=hold(i,j)
       enddo
       enddo
       go to 2000
      endif

c calculate the change of gradient
      absdg=0.d0
      do i=1,n3
       delg(i)=g(i)-gold(i)
       absdg=absdg+delg(i)**2
      enddo
      absdg=sqrt(absdg)

c calculate the error in the gradient
      do i=1,n3
       z(i)=delg(i)
       do j=1,n3
        z(i)=z(i)-hold(i,j)*dx(j)
       enddo
      enddo

c some dot products needed for the hessian perturbation
      dotzx=0.d0
      dotxx=0.d0
      dotdgx=0.d0
      dotzz=0.d0
      do i=1,n3
       dotzx=dotzx+z(i)*dx(i)
       dotxx=dotxx+dx(i)*dx(i)
       dotzz=dotzz+z(i)*z(i)
       dotdgx=dotdgx+delg(i)*dx(i)
      enddo
      absx=sqrt(dotxx)

c form three versions of the change in the hessian

c form hsr1
      do i=1,n3
      do j=1,n3
       hsr1(i,j)=z(i)*z(j)/dotzx
      enddo
      enddo
c form hpsb
      do i=1,n3
      do j=1,n3
       hpsb(i,j)=(dx(i)*z(j)+dx(j)*z(i))/dotxx -
     .            dotzx*dx(i)*dx(j)/dotxx**2      
      enddo
      enddo

c form hbfgs

      shs=0.d0
      do i=1,n3
      do j=1,n3
       shs=shs+dx(i)*hold(i,j)*dx(j)
      enddo
      enddo

      do i=1,n3
      do j=1,n3
       h2(i,j)=0.d0
       do k=1,n3
        h2(i,j)=h2(i,j)+hold(i,k)*hold(k,j)
       enddo
      enddo
      enddo

      do i=1,n3
      do j=1,n3
       hbfgs(i,j)=delg(i)*delg(j)/dotdgx - 
     .            dx(i)*h2(i,j)*dx(j)/shs
      enddo
      enddo

c      if(Nsearcxh.eq.3.or.Nsearch.eq.5)then
       if(negeigs.eq.0)then
c "Flow chart method" from Section 2.2.1 
c Theor Chem Acc (2016) 135:84

      ratio=dotzx/(absdg*absx)

      if(ratio.gt.0.1d0)then
       do i=1,n3
       do j=1,n3
        delh(i,j)=hbfgs(i,j)
       enddo
       enddo
      endif

      if(ratio.le.0.1d0)then
       do i=1,n3
       do j=1,n3
        delh(i,j)=hpsb(i,j)
       enddo
       enddo
      endif

      if(ratio.lt.-0.1d0)then
       do i=1,n3
       do j=1,n3
        delh(i,j)=hsr1(i,j)
       enddo
       enddo
      endif
c end the Nsearch if
      endif

c      if(nsearch.eq.4)then
      if(negeigs.gt.0)then
c TS approach of Bofill J Comp Chem 15, 1 (1994)
      phi=dotzx*dotzx/(dotzz*dotxx)
       do i=1,n3
       do j=1,n3
        delh(i,j)=phi*hsr1(i,j)+(1.d0-phi)*hpsb(i,j)
       enddo
       enddo
c end the Nsearch if
      endif

      do i=1,n3
      do j=1,n3
       h(i,j)=hold(i,j)+delh(i,j)
      enddo
      enddo
2000  continue

c now output revised combinedderivs
      open(unit=1,file='combinedderivs',status='old')
      jobtype=2
      write(1,50)com1
      write(1,*)jobtype
      write(1,50)com2
      write(1,*)natom
      write(1,50)com3
      do n=1,natom
      write(1,*)amas(n)
      enddo
      write(1,50)com4
      do n=1,natom
      write(1,51)lab(n)
      enddo
      write(1,50)com5
      write(1,*)energynew
      write(1,50)com6
      do n=1,natom
      write(1,33)(coord(k,n),k=1,3)
      enddo
      write(1,50)com7
      do n=1,natom
      do k=1,3
       nk=3*(n-1)+k
       write(1,33)g(nk)
      enddo
      enddo
      write(1,50)com8
        write(1,100)((h(i,j),j=1,i),i=1,n3)
      close(unit=1)

      end


      subroutine filelabel(k,ca)
      implicit double precision(a-h,o-z)

c this subroutines returns a character variable ca
c for the integer k

      character*20 ca
      character*1 I0(0:10),anum(10),ca1,ca2,ca3,ca4,ca5

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

      if(k.gt.9999.and.k.le.99999)then
      k1=k/10000
      k3=k/1000 - INT(k1)*10
      k4=k/100 - INT(k1)*100 - INT(k3)*10
      k5=k/10-INT(k1)*1000-k3*100-k4*10
      k6=k-INT(k1)*10000-k3*1000-k4*100-k5*10
      ca1=I0(k1)
      ca2=I0(k3)
      ca3=I0(k4)
      ca4=I0(k5)
      ca5=I0(k6)
      ca=ca1//ca2//ca3//ca4//ca5
      endif

      return
      end

