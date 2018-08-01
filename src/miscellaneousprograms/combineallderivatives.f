      program combineall

c this program combines the gradients (and hessians) from the several
c possible sources, into a single file containing the
c gradients and hessians

      implicit double precision(a-h,o-z)

      real*8, allocatable   :: g(:,:),h(:,:),g0(:,:),h0(:,:)

      real*8, allocatable   :: amas(:),bgh(:,:,:),coord(:,:)

      character*2, allocatable   :: lab(:)

      character*50 com1,com2,com3,com4,com5,com6,com7,com8


      open(unit=1,file='IN_CHARGES',status='old')
      read(1,*)
      read(1,*)ncharges
      close(unit=1)
c     ncharges > 0 if formal charges are present
c this means that some additional source of gradients will be present


c get the number of atoms
      open(unit=1,file='name.xyz',status='old')
      read(1,*)natom
      close(unit=1)

c allocate some necessary arrays
      allocate(g(natom,3))
c g is the combined gradient
      allocate(g0(natom,3))
c g0 is a temporary gradient
      allocate(h(3*natom,3*natom))
c h is the combined hessian
      allocate(h0(3*natom,3*natom))
c h0  is a temporary hessian
      allocate(lab(natom))
c lab is the elemental symbol for each atom
      allocate(amas(natom))
c amas is the atomic mass for each atom
      allocate(bgh(natom,3,3))
c bgh is used to temporarily store a hessian associated with
c embedded charges in fragment calculations
      allocate(coord(natom,3))
c coord are the Cartesin coordinates for each atom


      open(unit=1,file='combinedFRAGderivs',status='old')
c combinedFRAGderivs must be present
c combinedFRAGderivs contains the grdients (and hession) from the
c main fragment calculation at the input value of Level
      read(1,50)com1
      read(1,*)jobtype
c jobtype=0 energies only
c jobtype=1 gradients
c jobtype=2 hessians

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
      read(1,*)energy
      energyab=energy
c energy is the total energy from the main fragment calculations
      read(1,50)com6
      do n=1,natom
      read(1,33)(coord(n,k),k=1,3)
      enddo
33      format(1x,3D25.16)
      read(1,50)com7
      sumg=0.d0
      do n=1,natom
      read(1,33)g(n,1)
      read(1,33)g(n,2)
      read(1,33)g(n,3)
      sumg=sumg+g(n,1)**2+g(n,2)**2+g(n,3)**2
      enddo
c     write(6,*)'  Main Fragment gradient magnitude =      ',sqrt(sumg)
      read(1,50)com8
      n3=3*natom
100   format(6f10.6)
        read(1,100)((h(i,j),j=1,i),i=1,n3)

      close(unit=1)

c sym sec derivs
      do i=1,n3
      do j=1,i
       h(j,i)=h(i,j)
      enddo
      enddo

c add the ab initio nonbonded energy to the main fragment energy
c to get the total energy
      tote=0.d0
      open(unit=1,file='NearInt_energy',status='unknown')
       read(1,*,end=81)tote
81      close(unit=1)
      energy=energy+tote
      energyab=energyab+tote

c get the next set of derivatives
      open(unit=1,file='combinedABNBderivs',status='unknown')
c combinedABNBderivs may not be present if no ab initio nonbonded jobs
c were executed.
      read(1,*,end=30)
      read(1,*)
      sumg=0.d0
      do n=1,natom
       read(1,*)(g0(n,k),k=1,3)
       sumg=sumg+g0(n,1)**2+g0(n,2)**2+g0(n,3)**2
      enddo
c     write(6,*)' Ab initio nonbonded gradient magnitude = ',sqrt(sumg)
      if(jobtype.eq.2)then
        n3=3*natom
        read(1,100)((h0(i,j),j=1,i),i=1,n3)
      endif
      close(unit=1)
c sym sec derivs
      do i=1,n3
      do j=1,i
       h0(j,i)=h0(i,j)
      enddo
      enddo

c add to g and h
      do n=1,natom
      do k=1,3
       g(n,k)=g(n,k)+g0(n,k)
      enddo
      enddo
      do i=1,n3
      do j=1,n3
       h(i,j)=h(i,j)+h0(i,j)
      enddo
      enddo
30    continue


c get the next set
      open(unit=1,file='Electrostatic_derivs',status='unknown')
c Electrostatic_derivs contains derivatives from the long range
c electrostatics
      tote=0.d0
      read(1,*,end=32)tote

      energy=energy+tote
      energyelect=tote

      if(nflag.eq.0)go to 32

      read(1,*)
      read(1,*)
      sumg=0.d0
      do n=1,natom
       read(1,*)(g0(n,k),k=1,3)
       sumg=sumg+g0(n,1)**2+g0(n,2)**2+g0(n,3)**2
      enddo
c     write(6,*)' Electrostatic gradient magnitude       = ',sqrt(sumg)
      if(jobtype.eq.2)then
        n3=3*natom
        read(1,100)((h0(i,j),j=1,i),i=1,n3)
      endif
      close(unit=1)
c sym sec derivs
      do i=1,n3
      do j=1,i
       h0(j,i)=h0(i,j)
      enddo
      enddo
c add to g and h
      do n=1,natom
      do k=1,3
       g(n,k)=g(n,k)+g0(n,k)
      enddo
      enddo
      do i=1,n3
      do j=1,n3
       h(i,j)=h(i,j)+h0(i,j)
      enddo
      enddo
32    continue

      do n=1,natom
      do k=1,3
       write(10,*)g(n,k)
      enddo
      enddo

c get the next derivs
      open(unit=1,file='combinedFRAGbgderivs',status='unknown')
c combinedFRAGbgderivs contains derivatives due to embedded charges in
c the main fragment ab initio calculations
      read(1,*,end=90)
      read(1,*)
      sumg=0.d0
      do n=1,natom
       read(1,*)(g0(n,k),k=1,3)
       sumg=sumg+g0(n,1)**2+g0(n,2)**2+g0(n,3)**2
      enddo
c     write(6,*)' Fragmt embedded charges grad magnitude = ',sqrt(sumg)
      if(jobtype.eq.2)then
      do n=1,natom
       read(1,100)((bgh(n,k1,k2),k1=1,k2),k2=1,3)
      enddo
c sym sec derivs
      do k2=1,3
      do k1=1,k2
      do n=1,natom
       bgh(n,k2,k1)=bgh(n,k1,k2)
      enddo
      enddo
      enddo
      endif
      close(unit=1)

c add to g and h
      do n=1,natom
       do k=1,3
       g(n,k)=g(n,k)+g0(n,k)
      enddo
      enddo
      do n1=1,natom
      do k1=1,3
       j1=3*(n1-1)+k1
      do k2=1,3
       j2=3*(n1-1)+k2
       h(j1,j2)=h(j1,j2)+bgh(n1,k1,k2)
      enddo
      enddo
      enddo
90    continue


c get the next set

      open(unit=1,file='combinedABbgderivs',status='unknown')
c combinedABbgderivs contains derivatives due to embedded charges in
c ab initio nonbonded calculations
      read(1,*,end=31)
      read(1,*)
      sumg=0.d0
      do n=1,natom
       read(1,*)(g0(n,k),k=1,3)
       sumg=sumg+g0(n,1)**2+g0(n,2)**2+g0(n,3)**2
      enddo
c     write(6,*)' Nonbonded embed charges grad magnitude = ',sqrt(sumg)
      if(jobtype.eq.2)then
      do n=1,natom
       read(1,100)((bgh(n,k1,k2),k1=1,k2),k2=1,3)
      enddo
c sym sec derivs
      do k2=1,3
      do k1=1,k2
      do n=1,natom
       bgh(n,k2,k1)=bgh(n,k1,k2)
      enddo
      enddo
      enddo
      endif
      close(unit=1)
c add to g and h
      do n=1,natom
      do k=1,3
       g(n,k)=g(n,k)+g0(n,k)
      enddo
      enddo
      do n1=1,natom
      do k1=1,3
       j1=3*(n1-1)+k1
      do k2=1,3
       j2=3*(n1-1)+k2
       h(j1,j2)=h(j1,j2)+bgh(n1,k1,k2)
      enddo
      enddo
      enddo
31    continue

c get the double counted derivatives
      open(unit=1,file='Doublecount_derivs',status='unknown')
c Doublecount_derivs contains a correction to the gradient due
c to the twice countered interaction of embedded charges with other
c embedded charges

      tote=0.d0
      read(1,*,end=36)tote

      energy=energy+tote
      energyab=energyab+tote

      if(nflag.eq.0)go to 36

      read(1,*)
      read(1,*)
      sumg=0.d0
      do n=1,natom
       read(1,*)(g0(n,k),k=1,3)
       sumg=sumg+g0(n,1)**2+g0(n,2)**2+g0(n,3)**2
      enddo
c     write(6,*)' Doublecount gradient magnitude         = ',sqrt(sumg)
      if(jobtype.eq.2)then
        n3=3*natom
        read(1,100)((h0(i,j),j=1,i),i=1,n3)
      endif
      close(unit=1)
c sym sec derivs
      do i=1,n3
      do j=1,i
       h0(j,i)=h0(i,j)
      enddo
      enddo
c add to g and h
      do n=1,natom
      do k=1,3
       g(n,k)=g(n,k)+g0(n,k)
      enddo
      enddo
      do i=1,n3
      do j=1,n3
       h(i,j)=h(i,j)+h0(i,j)
      enddo
      enddo
36    continue


      do n=1,natom
      do k=1,3
       write(11,*)g(n,k)
      enddo
      enddo

c get the next set
      open(unit=1,file='Dispderivatives',status='unknown')
c Dispderivatives contains the derivative from the dispersion
c interaction
      dispenergy=0.d0
      read(1,*,end=70)dispenergy
      energy=energy+dispenergy

      if(nflag.eq.0)go to 70

      read(1,*)
      read(1,*)
      sumg=0.d0
      do n=1,natom
       read(1,*)(g0(n,k),k=1,3)
       sumg=sumg+g0(n,1)**2+g0(n,2)**2+g0(n,3)**2
      enddo
c     write(6,*)' Dispersion gradient magnitude          = ',sqrt(sumg)
      if(jobtype.eq.2)then
        n3=3*natom
        read(1,100)((h0(i,j),j=1,i),i=1,n3)
      endif
      close(unit=1)
c sym sec derivs
      do i=1,n3
      do j=1,i
       h0(j,i)=h0(i,j)
      enddo
      enddo

c add to g and h
      do n=1,natom
      do k=1,3
       g(n,k)=g(n,k)+g0(n,k)
      enddo
      enddo
      do i=1,n3
      do j=1,n3
       h(i,j)=h(i,j)+h0(i,j)
      enddo
      enddo

70    continue

c get the next set
      open(unit=1,file='Induction_derivs',status='unknown')
c Induction_derivs contains the derivative from the 
c induction calculation
      tote=0.d0
      read(1,*,end=34)Eind

      energy=energy+Eind

      if(nflag.eq.0)go to 34

      read(1,*)
      read(1,*)
      sumg=0.d0
      do n=1,natom
       read(1,*)(g0(n,k),k=1,3)
       sumg=sumg+g0(n,1)**2+g0(n,2)**2+g0(n,3)**2
      enddo
c     write(6,*)' Induction gradient magnitude           = ',sqrt(sumg)
      if(jobtype.eq.2)then
        n3=3*natom
        read(1,100)((h0(i,j),j=1,i),i=1,n3)
      endif
      close(unit=1)
c sym sec derivs
      do i=1,n3
      do j=1,i
       h0(j,i)=h0(i,j)
      enddo
      enddo
c add to g and h
      do n=1,natom
      do k=1,3
       g(n,k)=g(n,k)+g0(n,k)
      enddo
      enddo
      do i=1,n3
      do j=1,n3
       h(i,j)=h(i,j)+h0(i,j)
      enddo
      enddo
34    continue

      go to 1000



c get the next set


c     go to 80

      open(unit=1,file='Crosschargederivatives',status='unknown')
      read(1,*,end=80)
      read(1,*)
      sumg=0.d0
      do n=1,natom
       read(1,*)(g0(n,k),k=1,3)
       sumg=sumg+g0(n,1)**2+g0(n,2)**2+g0(n,3)**2
      enddo
c     write(6,*)' Crosscharge gradient = ',sqrt(sumg)
      if(jobtype.eq.2)then
        n3=3*natom
        read(1,100)((h0(i,j),j=1,i),i=1,n3)
      endif
      close(unit=1)
c sym sec derivs
      do i=1,n3
      do j=1,i
       h0(j,i)=h0(i,j)
      enddo
      enddo

c add to g and h
      do n=1,natom
      do k=1,3
       g(n,k)=g(n,k)+g0(n,k)
      enddo
      enddo
      do i=1,n3
      do j=1,n3
       h(i,j)=h(i,j)+h0(i,j)
      enddo
      enddo

80    continue
c finally get the ab initio nb energy
c      open(unit=1,file='NearInt_energy',status='unknown')
c      abenergy=0.d0
c      read(1,*,end=60)abenergy
c60    continue
c      close(unit=1)

c      energy=energy+abenergy

1000  continue
c now all the derivatives are combined

c output energies to SMFA.out
      open(unit=1,file='SMFA.out.temp',status='unknown')
      write(1,*)' *****************************************************'
      write(1,*)
      write(1,*)' The total electronic energy of the molecule is'
      write(1,*)
      write(1,600)energy
600   format(f25.8)
      write(1,*)
      write(1,*)' *****************************************************'
      write(1,*)' This energy is composed of an ab initio component of'
      write(1,600)energyab
      if(abs(energyelect).gt.1.d-8)then
      Eelect=energyelect*2625.5d0
      write(1,*)' and a long range electrostatic component of'
      write(1,601)energyelect,' = ',Eelect,' (kJ/mol)'
      endif
601   format(f25.8,a3,f12.4,a9)
      if(abs(Eind).gt.1.d-9)then
      Einduction=Eind*2625.5d0
      write(1,*)' and an induction energy of'
      write(1,601)Eind,' = ',Einduction,' (kJ/mol)'
      endif
      if(abs(dispenergy).gt.1.d-8)then
      Edisp=dispenergy*2625.5d0
      write(1,*)' and a long range dispersion energy of'
      write(1,601)dispenergy,' = ',Edisp,' (kJ/mol)'
      endif
      close(unit=1)


      sumg=0.d0
      do n=1,natom
       sumg=sumg+g(n,1)**2+g(n,2)**2+g(n,3)**2
      enddo
c     write(6,*)
c     write(6,*)' Total gradient magnitude               = ',sqrt(sumg)

c we ouput them in the same format as the first file

      open(unit=1,file='combinedderivs',status='unknown')
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
      write(1,*)energy
      write(1,50)com6
      do n=1,natom
      write(1,33)(coord(n,k),k=1,3)
      enddo
      write(1,50)com7
      do n=1,natom
      write(1,33)g(n,1)
      write(1,33)g(n,2)
      write(1,33)g(n,3)
      enddo
      write(1,50)com8
      n3=3*natom
        write(1,100)((h(i,j),j=1,i),i=1,n3)

      do n=1,natom
      write(12,33)g(n,1)
      write(12,33)g(n,2)
      write(12,33)g(n,3)
      enddo

      close(unit=1)
   
      end 
