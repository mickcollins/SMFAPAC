      program optimiser   
      implicit double precision(a-h,o-z)

      real*8, allocatable   :: g(:,:),h(:,:),eig(:),v(:,:),hs(:,:)
      real*8, allocatable   :: gc(:,:),hc(:,:)
      real*8, allocatable   :: coeff(:),dis(:,:),ds(:)

      real*8, allocatable   :: amas(:),coord(:,:),work(:),dot(:),x(:,:)

      character*2, allocatable   :: lab(:)

      character*50 com1,com2,com3,com4,com5,com6,com7,com8

      character*1 jobz, uplo
      character*3 fl1,fl2,fl3,fl4
C
C-------------------------------------------------------------------------------
C    The input is the same as Mick's code except that some of the arrays
C       are stored differently
C    This does not affect the files read in, nor those written at the end
C------------------------------------------------------------------------------

      Bohr=1.88972598860

      open(unit=1,file='combinedderivs',status='old')
      read(1,50)com1
      read(1,*)jobtype
      read(1,50)com2
      read(1,*)natom

      allocate(g(3,natom))
      allocate(gc(3,natom))
      n3=3*natom
      allocate(h(n3,n3))
      allocate(hc(n3,n3))
      allocate(hs(n3,n3))
      allocate(eig(n3))
      allocate(coeff(n3))
      allocate(dot(n3))
      allocate(v(n3,n3))
      allocate(amas(natom))
      allocate(coord(3,natom))
      allocate(x(3,natom))
      allocate(dis(3,natom))
      allocate(ds(natom))
      allocate(lab(natom))
      allocate(work(3*n3))

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
      read(1,*)energy
      write(6,*)' The input energy = ',energy
      read(1,50)com6
      do n=1,natom
      read(1,33)(coord(k,n),k=1,3)
      enddo
c coord is in Bohr
33      format(1x,3D25.16)
      read(1,50)com7
      do n=1,natom
      read(1,*)g(1,n)
      read(1,*)g(2,n)
      read(1,*)g(3,n)
      enddo
      read(1,50)com8
      n3=3*natom
100   format(6f10.6)
        read(1,100)((h(i,j),j=1,i),i=1,n3)

      close(unit=1)

c the gradient is in Hartree/ Bohr
c the hessian is in Hartree / Bohr**2
c the coordinates are in Bohr

c sym
      do i=1,n3
      do j=1,i
      h(j,i)=h(i,j)
      enddo
      enddo

c store h
      do i=1,n3
      do j=1,n3
      hs(i,j)=h(i,j)
      enddo
      enddo
c find the number of negative eigenvalues, as this is needed
c by the hessian update program
      call  project(hs,coord,natom)

      lwork=3*n3
      call dsyev('N','U',n3,hs,n3,eig,work,lwork,ifail)
      nneg=0
      do i=1,n3
       if(eig(i).lt.-1.d-4)nneg=nneg+1
      enddo
      open(unit=8,file='NEGEIGVALS',status='unknown')
      write(8,*)nneg
      close(unit=8)


c read in the stepnumber in the optiisation
      open(unit=2,file='STEPNUMBER',status='old')
       read(2,*)istep
        write(6,*)' Optimisation step number = ',istep
      close(unit=2)

c read in whether to search for a minimum or a TS
      open(unit=2,file='IN_JOBTYPE_original',status='old')
      read(2,*)
      read(2,*)Nsearch
      close(unit=2)
      if(Nsearch.eq.4)write(6,*)' TS search'

c get the gradients and hessians for any constraints that are input
c in IN_CONSTRAINTS
c     gc=0.d0
c     hc=0.d0
c     call optconstraints(natom,coord,value,gc,hc,jobtype)

c add on the constraint foreces and hession
c     do n=1,natom
c     do k=1,3
c      write(34,*)g(n,k),gc(n,k),abs(gc(n,k)-g(n,k))
c      g(n,k)=g(n,k)+gc(n,k)
c     enddo
c     enddo
c     do i=1,n3
c     do j=1,n3
c   use hessian of the constraints
c      h(i,j)=h(i,j)+hc(i,j)
c     enddo
c     enddo

c      project translations and rotations out of hessian

      call  project(h,coord,natom)

c      Trust radius
c      to get this to work need to store previous energy and
c      previous trust radius
c      set to temporary value for now

       open(unit=8,file='TRUSTRADIUS',status='old')
       read(8,*)radius
       close(unit=8)
    
      call constraints(natom,coord,cvalue,gc,hc,jobtype) 

      do k=1,3
      do n=1,natom
       g(k,n)=g(k,n)+gc(k,n)
      enddo
      enddo

      do i=1,n3
      do j=1,n3
       h(i,j)=h(i,j)+hc(i,j)
      enddo
      enddo

      if(Nsearch.eq.4)then
      write(67,*)' Calling hesaugTS'
       call hesaugTS(n3,h,g,dis,radius,energy,istep)
      else
       call hesaug(n3,h,g,dis,radius,energy,istep)
      endif

c estimate the change in energy (just the linear change)

c first remove the constraint parts

      do k=1,3
      do n=1,natom
       g(k,n)=g(k,n)-gc(k,n)
      enddo
      enddo

      do i=1,n3
      do j=1,n3
       h(i,j)=h(i,j)-hc(i,j)
      enddo
      enddo


      deltae=ddot(n3,dis,1,g,1)


c     dstep=ddot(n3,dis,1,dis,1)
      

c add the quadratic part
      esum=0.d0
      do k1=1,3
      do k2=1,3
      do n1=1,natom
      do n2=1,natom
       i=3*(n1-1)+k1
       j=3*(n2-1)+k2
       esum=esum+dis(k1,n1)*h(i,j)*dis(k2,n2)
      enddo
      enddo
      enddo
      enddo

      open(unit=8,file='predictedenergychange',status='unknown')
       write(8,*)deltae+0.5d0*esum
       totaldeltaE=deltae+0.5d0*esum
      close(unit=8)

c     update coordinates 

c if istep = 0, leave the coordinates unchanged, as the gradient is
c only a SCF gradient

      if(istep.eq.0)then
      do n=1,natom
      do k=1,3
       dis(k,n)=0.d0
      enddo
      enddo
      endif

      do n=1,natom
      do k=1,3
       x(k,n)=(coord(k,n)+dis(k,n))/bohr
      enddo
      enddo

c output the change in coordinates in Bohr
      open(unit=3,file='deltax',status='unknown')
      do n=1,natom
      do k=1,3
       write(3,*)dis(k,n)
      enddo
      enddo
      close(unit=3)

c if istep = 0, skip to the end
      if(istep.eq.0)go to 500

c ---------------------------------------------------
c         the tests beyond here look fairly arbitrary
C         these are the same as Mick's code except that
c         they have all been moved to the end to be together
c---------------------------------------------------

c estimate the step in coordinates
      step=0.d0
      do n=1,natom
      do k=1,3
       step=step+dis(k,n)**2
      enddo
      enddo
      step=sqrt(step)/Bohr

c add the constraint forces back in
      do k=1,3
      do n=1,natom
       g(k,n)=g(k,n)+gc(k,n)
      enddo
      enddo



c calculate thresold quantities
      rmsg=0.d0
      do n=1,natom
      do k=1,3
       rmsg=rmsg+g(k,n)**2
      enddo
      enddo
       rmsg=sqrt(rmsg/n3)

c calc max gradient
      gmax=0.d0
      do n=1,natom
      do k=1,3
       if(abs(g(k,n)).gt.gmax)gmax=abs(g(k,n))
      enddo
      enddo


      rmsd=0.d0
      stepmax=0.d0
      do n=1,natom
      do k=1,3
      if(abs(dis(k,n)).gt.stepmax)stepmax=abs(dis(k,n))
      rmsd=rmsd+dis(k,n)**2
      enddo
      enddo
      rmsd=sqrt(rmsd/n3)


      write(6,*)' The current gradient norm = ',rmsg 
c     write(6,*)' The accessible gradient = ',goptsize
      write(6,*)' The estimated next step size = ',step
c     write(6,*)' The estimated change in energy would be = ',deltae
      write(6,*)' The estimated change in energy  = ',totaldeltaE

c     these tolerences may need adjusting

      rmsgtol=1.5d-4
      gmaxtol=4.5d-4
      rmsdtol=1.2d-3
      stepmaxtol=1.8d-3

      if(gmax.lt.gmaxtol)then
       fl1='YES'
      else
       fl1='NO'
      endif
      if(rmsg.lt.rmsgtol)then
       fl2='YES'
      else
       fl2='NO'
      endif
      if(stepmax.lt.stepmaxtol)then
       fl3='YES'
      else
       fl3='NO'
      endif
      if(rmsd.lt.rmsdtol)then
       fl4='YES'
      else
       fl4='NO'
      endif

      write(6,*)'          Item               Value      Thresold'
      write(6,150)'Maximum Force       ',gmax,gmaxtol,fl1
      write(6,150)'RMS     Force       ',rmsg,rmsgtol,fl2
      write(6,150)'Maximum Displacement',stepmax,stepmaxtol,fl3
      write(6,150)'RMS     Displacement',rmsd,rmsdtol,fl4
150   format(a20,5x,2f12.6,3x,a3)


c look for convergence

      nscore=0
      if(fl1.eq.'YES')nscore=nscore+1
      if(fl2.eq.'YES')nscore=nscore+1
      if(fl3.eq.'YES')nscore=nscore+1
      if(fl4.eq.'YES')nscore=nscore+1

      nscoref=0
      if(fl1.eq.'YES')nscoref=nscoref+1
      if(fl2.eq.'YES')nscoref=nscoref+1

500   continue

c output the new coordinates
      open(unit=2,file='Newcoords.xyz',status='unknown')
      write(2,*)natom
      write(2,*)' coordinates output from Opt'
      do n=1,natom
       write(2,200)lab(n),(x(k,n),k=1,3)
       write(20,200)lab(n),(x(k,n),k=1,3)
      enddo
c     write(2,*)'END'
      close(unit=2)      

200   format(a2,3f15.8)


      if(nscore.ge.3.and.istep.gt.0)then
       write(6,*)' Geometry has converged'
       stop
      endif
      
c     if(nscoref.ge.2)then
c      write(6,*)' Geometry has converged'
c      stop
c     endif

      end

      subroutine project(fcm,c,nat)
      implicit double precision (a-h,o-z)

C     uses the Miller,Handy and Adams algorithm
C     to project out rotational and translational
C     components of teh force constant matrix
C     (strickly only valid at a stationary point
C      but can be useful in geometry optimisations)

C     fcm = force constant matrix
C     C = cartesian geometry
C     Nat = number of atoms

C     Caution - overwrites original force constant matrix
C     fcm with projected version

      dimension c(3,nat),fcm(nat*3,nat*3)
c     temporary arrays
      dimension f1(nat*3,nat*3),al(nat*3,nat*3)

      nat3=nat*3
      call zero(al,nat3*6)
      DO  N=1,NAT
      AL(3*N-2,1)=1.0D0
      AL(3*N-1,2)=1.0D0
      AL(3*N,3)=1.0D0
      AL(3*N-2,4)=-C(2,N)
      AL(3*N-1,4)=C(1,N)
      AL(3*N-1,5)=-C(3,N)
      AL(3*N,5)=C(2,N)
      AL(3*N-2,6)=C(3,N)
      AL(3*N,6)=-C(1,N)
      ENDDO

      call  orthy(al,nat3,6)
      call dgemm('N','T',nat3,nat3,6,-1d0,al,nat3,al,nat3,0d0,f1,nat3)

      do i=1,nat3
      f1(i,i)=f1(i,i)+1d0
      enddo

C
      call dgemm('N','N',nat3,nat3,nat3,1d0,f1,nat3,fcm,nat3,0d0,al,
     1 nat3)
      call dgemm('N','N',nat3,nat3,nat3,1d0,al,nat3,f1,nat3,0d0,fcm,
     1   nat3)

      RETURN
      END
      subroutine orthy(a,m,n)
      implicit double precision (a-h,o-z)
      dimension a(m,n), work(n)

c     orthonormalises the columns of matrix A
c     caution : intended as part of subroutine projec
c     not suitable for general use
c     sets small elements of A to zero
c     small aribitrarily defined as 1d-10

      parameter (small=1d-10)

      do j=1,n
      do i=1,m
      if(dabs(a(i,j)).lt.small)a(i,j)=0d0
      enddo
      enddo

      c=dsqrt(ddot(m,a(1,1),1,a(1,1),1))

      if(c.gt.small)then
       a(1:m,1)=a(1:m,1)/c
      else
       a(1:m,1)=0d0
      endif


      do na=2,n

          do nb=1,na-1
          work(nb)=ddot(m,a(1,na),1,a(1,nb),1)
          enddo

          do nb=1,na-1
          call daxpy(m,-work(nb),a(1,nb),1,a(1,na),1)
          enddo

          do i=1,m
          if(dabs(a(i,na)).lt.small)a(i,na)=0d0
          enddo

          c=dsqrt(ddot(m,a(1,na),1,a(1,na),1))
          if(c.gt.small)then
          a(1:m,na)=a(1:m,na)/c
          else
          a(1:m,na)=0d0
          endif

      enddo
      return
      end
      subroutine zero(a,n)
      double precision a(n)
      do i=1,n
      a(i)=0d0
      enddo
      return
      end
      subroutine hesaug(n,fcm,grad1,dgeom1,radius,energy,istep)

      implicit double precision (a-h,o-z)
      dimension grad1(3,n/3),dgeom1(3,n/3)
      dimension fcm(n,n),grad(n),dgeom(n)
      dimension a(n+1,n+1),v(n+1,n+1),w(n+1),w2(n+1)
      dimension work(10*n)

      character*20 ca,ca1

      parameter (small=2d-5)

      do k=1,3
      do i=1,n/3
       j=3*(i-1)+k
       grad(j)=grad1(k,i)
      enddo
      enddo

      lwork=10*n

c     there will be 6 eigenvalues very close to zero since the 
c       translations and rotations are projected out
c      however there may be some noise and do not want to
c      search along very flat direction, so test for eigenvalues
c      but DO want to search along hydrogen bonds so this
c      threshold is set to 0.0001 which is approx 25 wavenumbers    

c if appropriate get the last predicted energy change
      if(istep.gt.0)then
       open(unit=8,file='predictedenergychange',status='old')
       read(8,*)predicted
       close(unit=8)
       laststep=istep-1
       call filelabel(laststep,ca1)
       n1=index(ca1,' ')-1
       ca='energy.'//ca1(1:n1)
       open(unit=8,file=ca,status='old')
       read(8,*)oldenergy
       close(unit=8)

       actualdeltaE=energy-oldenergy
       write(6,*)' Previously predicted energy change = ',predicted
       write(6,*)' Actual energy change               = ',actualdeltaE
       write(6,*)' The trust radius now = ',radius
      endif


      do j=1,n
      do i=1,n
      v(i,j)=fcm(i,j)
      enddo
      enddo

      call dsyev('V','U',n,v,n+1,w,work,lwork,ifail)

c     write(6,*) 'eigenvalues of projected hessian'
c     write(6,10) (w(i),i=1,n)
10    format(5f15.5)


c     count the number of non-zero values

      nz=0
      nzero=0
      do i=1,n
      if(dabs(w(i)).lt.small)then
      w(i)=0d0
      nzero=nzero+1
      endif
      enddo
c     write(6,*) 'There are ', nzero  , ' zero eigenvalues'
      nz=n-nzero

      nneg=0
      do i=1,n
      if(w(i).lt.-small) nneg=nneg+1
      enddo
      write(6,*) 'The hessian has ', nneg ,' negative eigenvalues'


c    though we have eliminates very flat eigenvectors, 
c    do not want to take large steps so shift all eigenvalues
c    up, so lowest is 0.05
      shift=0.05d0
      if(nneg.ne.0) then
      shift=-w(1)+0.05d0
      endif
c     write(6,*) 'apply level shift ',shift

      do j=1,n+1
      do i=1,n+1
      a(i,j)=0d0
      enddo
      enddo

c     set up augmented hessian in normal mode basis
c     put eigenvalues of fcm along diagonal
c     ignoring those which are zero
c     and applying a level shift so that negative
c     eigenvalues become positive
c     NOTE : this assumes are looking for a minimum
c     in the coordinate space. This could be altered
c     to look for other types of stationary points
c see Anglada and Bofill, Int J Quant Chem 62,153 (1997)

      m=0
      do i=1,n
      if(w(i).ne.0d0)then
      m=m+1
      a(m,m)=w(i)+shift
      a(nz+1,m)=ddot(n,v(1,i),1,grad,1)
      a(m,nz+1)=a(nz+1,m)
      endif
      enddo

c     cannot remember why this line was there
c yes, should be here, see
c Anglada and Bofill, Int J Quant Chem 62,153 (1997) Eq 12
      a(nz+1,nz+1)=-shift*radius*radius


      call dsyev('V','U',nz+1,a,n+1,w2,work,lwork,ifail)

c     write(6,*)'eigenvalues of augmented hessian'
c     write(6,10) (w2(i),i=1,nz+1)

c     calculate the displacement

      do i=1,nz
      work(i)=a(i,1)/a(nz+1,1)
      enddo
      do i=1,n
      dgeom(i)=0d0
      enddo
      m=0
      do i=1,n
      if(w(i).ne.0d0)then
      m=m+1
      do j=1,n
      dgeom(j)=dgeom(j)+work(m)*v(j,i)
      enddo
      endif
      enddo

      gnorm=dsqrt(ddot(n,dgeom,1,dgeom,1))

c     scale using trust radius
c     NB the trust radius should really be set dynamically
c simpler than Anglada and Bofill
      if(istep.gt.1)then
       ratio=actualdeltaE/predicted
       if(gnorm.gt.0.95d0*radius)then
        if(ratio.gt.0.85.and.ratio.lt.1.15)then
         radius=radius*1.4142d0
        else
         if(ratio.gt.0.75.and.ratio.lt.1.25)radius=radius/2.d0
        endif
       endif

       open(unit=8,file='TRUSTRADIUS',status='old')
        write(8,*)radius
       close(unit=8)
      endif

      if(gnorm.gt.radius)then
      do i=1,n
      dgeom(i)=dgeom(i)*radius/gnorm
      enddo
      endif

      do k=1,3
      do i=1,n/3
       j=3*(i-1)+k
       dgeom1(k,i)=dgeom(j)
      enddo
      enddo

      return
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

