      program frequencies
      implicit double precision(a-h,o-z)

      real*8, allocatable   :: g(:,:),h(:,:),eig(:),v(:,:)

      real*8, allocatable   :: amas(:),coord(:,:),work(:),freqw(:)

      real*8, allocatable   :: dipdr(:,:,:),height(:),weight(:,:)

      character*2, allocatable   :: lab(:)

      character*50 com1,com2,com3,com4,com5,com6,com7,com8

      character*1 jobz, uplo

      character*2 ele

      open(unit=1,file='combinedderivs',status='old')
      read(1,50)com1
      read(1,*)jobtype
      read(1,50)com2
      read(1,*)natom

      allocate(g(3,natom))
      n3=3*natom
      allocate(h(n3,n3))
      allocate(eig(n3))
      allocate(v(n3,n3))
      allocate(amas(natom))
      allocate(coord(3,natom))
      allocate(lab(natom))
      allocate(work(9*natom))
      allocate(dipdr(3,3,natom))
      allocate(weight(3,n3))
      allocate(height(n3))
      allocate(freqw(n3))

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
      read(1,50)com6
      do n=1,natom
      read(1,33)(coord(k,n),k=1,3)
      enddo
33      format(1x,3D25.16)
      read(1,50)com7
      do n=1,natom
      read(1,33)g(1,n)
      read(1,33)g(2,n)
      read(1,33)g(3,n)
      enddo
      read(1,50)com8
      n3=3*natom
100   format(6f10.6)
        read(1,100)((h(i,j),j=1,i),i=1,n3)

      close(unit=1)


c sym
      do i=1,n3
      do j=1,i
       h(j,i)=h(i,j)
      enddo
      enddo

c the default width (FWHM) is 5.0 cm-1
      width=5.d0

c we look to see if IN_ISOTOPES is present
c If it is, then we change the masses in one of two ways
      open(unit=1,file='IN_ISOTOPES',status='unknown')
      read(1,*,end=500)nflag
c here because isotopes are active
c if nflag=1, all atoms of a given element are changed
c if nflag=2, we read atom numbers to change
      if(nflag.eq.1)then
       read(1,*)ele,tmass
       do n=1,natom
        if(lab(n).eq.ele)amas(n)=tmass
       enddo
       read(1,*)width
       go to 500
      endif
      if(nflag.eq.2)then
       read(1,*)nchanges
       if(nchanges.gt.0)then
       do i=1,nchanges
        read(1,*)n1,am1
        amas(n1)=am1
       enddo
       endif
       read(1,*)width
      endif 

500    close(unit=1)



c project out rotation and translation
      call project(h,coord,natom)

c     do i=1,n3
c     do j=1,i
c      write(55,*)h(i,j)
c     enddo
c     enddo

c mass weight
      do n1=1,natom
      do n2=1,natom
      do k1=1,3
      do k2=1,3
       i1=3*(n1-1)+k1
       i2=3*(n2-1)+k2
       h(i1,i2)=h(i1,i2)/sqrt(amas(n1)*amas(n2))
      enddo
      enddo
      enddo
      enddo

c eigs only
c     jobz='N'
c eigs and vectors
      jobz='V'
      uplo='U'
      lwork=3*n3

      call dsyev(jobz,uplo,n3,h,n3,eig,work,lwork,info)
      if(info.ne.0)then
       write(6,*)' error in dsyev ',info
       stop
      endif

c h now contains the eigenvectors

c read in the dipole moment derivatives
      open(unit=1,file='combDipDerivs',status='old')
      do n=1,natom
      do k=1,3
       read(1,*)(dipdr(j,k,n),j=1,3)
      enddo
      enddo  
      close(unit=1)
    
c mass weight the derivatives
      do n=1,natom
      do k=1,3
      do j=1,3
       dipdr(j,k,n)=dipdr(j,k,n)/sqrt(amas(n))
      enddo
      enddo
      enddo

c calculate the magnitude of the dipole derivative wrt the normal coord
      do j=1,3
      do i=1,n3
       weight(j,i)=0.d0
       do n=1,natom
       do k=1,3
        nk=3*(n-1)+k
        weight(j,i)=weight(j,i)+dipdr(j,k,n)*h(nk,i)
       enddo
       enddo
      enddo
      enddo

      hmax=0.d0
      do i=1,n3
       height(i)=0.d0
       do j=1,3
        height(i)=height(i)+weight(j,i)**2
       enddo
       if(height(i).gt.hmax)hmax=height(i)
      enddo

c     do i=1,n3
c      height(i)=height(i)/hmax
c     enddo

      H2wav=5140.5

      open(unit=1,file='FREQUENCIES',status='unknown')
      write(1,*)" Frequency (cm-1)       Intensity (arb. units)"
      zpe=0.d0
      do i=1,n3
      if(abs(eig(i)).lt.1.d-12)eig(i)=1.d-12
      esign=eig(i)/abs(eig(i))
      freqw(i)=esign*sqrt(abs(eig(i)))*H2wav
      if(eig(i).gt.1.d-11)zpe=zpe+sqrt(eig(i))*61.494d0
c      write(1,200)freqw(i),height(i)
       write(1,200)freqw(i),height(i)*freqw(i)
      enddo
      write(1,*)
      write(1,*)
      write(1,222)' The zero-point energy = ',0.5d0*zpe,' kJ/mol'
222   format(a25,f10.2,a7)
      close(unit=1)
      open(unit=1,file='NORMAL_MODES',status='unknown')
      write(1,*) '       Cartesian displacements for each normal mode'
      do i=1,n3
       write(1,*)" Mode ",i,"   Frequency (cm-1) ",freqw(i)
       write(1,*)
       write(1,*)"   Atom            x           y          z"
       do n=1,natom
        n1=3*(n-1)
        write(1,101)n,(h(n1+j,i)/sqrt(amas(n)),j=1,3)
       enddo
      enddo
      close(unit=1)
101   format(i6,6x,3f12.4)
c calculate and write out the spectrum
      open(unit=1,file='SPECTRUM',status='unknown')
      write(1,*)" Frequency (cm-1)    Intensity (arb. units)"
      write(1,*)

      width=width/2.d0
      f=100.d0
      do m=1,9800
       f=f+0.5d0
       psum=0.d0
       do i=1,n3
        x=(f-freqw(i))/width
c       psum=psum+height(i)/(1.d0+x**2)
        psum=psum+freqw(i)*height(i)/(1.d0+x**2)
       enddo
       write(1,102)f,psum
      enddo
      close(unit=1)

102   format(f12.1,10x,f12.6)


200   format(f12.1,12x,f12.4)



      end


