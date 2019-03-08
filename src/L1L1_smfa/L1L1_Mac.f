
      program L1L1_mac

      implicit double precision(a-h,o-z)

      parameter(nsmall=60)

c nsmall is the same parameter in SMFA, the maximum number of
c groups in a fragment.

      integer, allocatable   ::ksign(:),numat(:),natstore(:,:)

      integer, allocatable   :: natbig(:),numbig(:,:)

      integer, allocatable   ::matstore(:,:),nsign(:)

      integer, allocatable   ::numchggrps(:),chchggrps(:)

      integer, allocatable   ::nom(:,:)

      integer, allocatable   :: n4ab(:)

      dimension mata(2,2)

c     dimension rad(91)
c     dimension numb(91)
      dimension rad(110), am(110)
      dimension numb(110)

      integer, allocatable  ::nat0(:),nch0(:)
      real*8, allocatable   ::c0(:,:,:),radius(:,:),dmin(:,:)
      integer, allocatable  ::nat1(:),nch1(:)
      real*8, allocatable   ::c1(:,:,:)
      character*2, allocatable  ::lab0(:,:),lab1(:,:)

      character*80 abinitio(6),inputformat

      character*20 ca,ca1,ca2

c     character*2 sym(91)
      character*2 sym(110)

      data (sym(i),i=1,110)/'H','He','Li','Be','B','C','N','O','F','Ne',
     .'Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr',
     .'Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb',
     . 'Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',
     . 'Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu',
     . 'Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os',
     . 'Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac',
     . 'Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No',
     . 'Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds'/
      data (am(i),i=1,110)/1.00794,4.002602,6.941,9.012182,10.811,
     $     12.0107,14.0067,15.9994,18.9984032,20.1797,22.98977,24.305,
     $     26.981538,28.0855,30.973761,32.065,35.453,39.948,39.0983,
     $     40.078,44.95591,47.867,50.9415,51.9961,54.938049,55.845,
     $     58.9332,58.6934,63.546,65.409,69.723,72.64,74.9216,78.96,
     $     79.904,83.798,85.4678,87.62,88.90585,91.224,92.90638,95.94,
     $     98,101.07,102.9055,106.42,107.8682,112.411,114.818,118.71,
     $     121.76,127.6,126.90447,131.293,132.90545,137.327,138.9055,
     $     140.116,140.90765,144.24,145,150.36,151.964,157.25,158.92534,
     $     162.5,164.93032,167.259,168.93421,173.04,174.967,178.49,
     $     180.9479,183.84,186.207,190.23,192.217,195.078,196.96655,
     $     200.59,204.3833,207.2,208.98038,209,210,222,223,226,227,
     $     232.0381,231.03588,238.02891,237,244,243,247,247,251,252,257,
     $     258,259,262,261,262,266,264,277,268,281/
      data (rad(i),i=1,110)/
     .  0.23, 1.22, 0.68, 0.35, 0.83, 0.68, 0.68, 0.68, 0.64, 1.6,
     .  0.97, 1.1,  1.35, 1.2,  1.05, 1.02, 0.99, 1.92, 1.33, 0.99,
     .  1.44, 1.47, 1.33, 0.67, 1.35, 1.34, 1.33, 1.5,  1.52, 1.45,
     .  1.22, 1.17, 1.21, 1.22, 1.21, 1.98, 1.47, 1.12, 1.78, 1.56,
     .  1.48, 1.47, 1.35, 1.4,  1.45, 1.5,  1.59, 1.69, 1.63, 1.46,
     .  1.46, 1.47, 1.4,  2.18, 1.67, 1.34, 1.87, 1.83, 1.82, 1.81,
     .  1.8,  1.8,  1.99, 1.79, 1.76, 1.75, 1.74, 1.73, 1.72, 1.94,
     .  1.72, 1.57, 1.43, 1.37, 1.35, 1.37, 1.32, 1.5,  1.5,  1.7,
     .  1.55, 1.54, 1.54, 1.68, 1.21, 1.50, 1.50, 1.9,  1.88, 1.79,
     .  1.61, 1.58, 1.55, 1.53, 1.51, 0.99, 1.54, 1.83, 1.50, 1.50,
     .  1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50, 1.50/


      bigtol=15.d0
      bigtol2=100.d0
100   format(a2,3f13.6)

c      open(unit=7,file='atomic.dat',status='old')
c      do m=1,91
c      read(7,*)numb(m)
c      read(7,212)sym(m)
c      read(7,*)amm
c      read(7,*)rad(m)
c      enddo
c212   format(a2)
c      close(unit=7)

c read in a parameter that allows or prevents caps
        open(unit=1,file='IN_CAPS',status='old')
        read(1,*)
        read(1,*)nocapsatall
c c if nocapsatall=1, there are no caps at all
        close(unit=1)

      do m=1,110
       numb(m)=m
      enddo

90     format(a80)
c read in some job control parameters
c need dtol
      open(unit=1,file='IN_FRACT',status='old')
      read(1,90)icomm
      read(1,*)Level
      read(1,90)icomm
      read(1,*)dtol
      close(unit=1)

      open(unit=1,file='frags.out_Lev0',status='old')
      read(1,*)
      read(1,*)ngroups
      allocate(nat0(ngroups))
      allocate(nch0(ngroups))
      allocate(c0(ngroups,60,3))
      allocate(lab0(ngroups,60))
      allocate(radius(ngroups,60))
c     allocate(numa(ngroups,60))
      allocate(dmin(ngroups,ngroups))
      allocate(nom(ngroups,ngroups))

      close(unit=1)

c output the number of groups, as the script "Run_polar.pl" needs it
      open(unit=8,file='OUT_NGROUPS',status='unknown')
       write(8,*)ngroups
      close(unit=8)

c read in the Level 0 (group) coordinates
      do n=1,ngroups
       call filelabel(n,ca)
       ca1='Lev0_COORD'//ca
       open(unit=1,file=ca1,status='old')
       read(1,*)nch0(n)
       m=1
1      continue
c      read(1,100,end=2)lab0(n,m),(c0(n,m,k),k=1,3)
       read(1,*,end=2)lab0(n,m),(c0(n,m,k),k=1,3)
       m=m+1
       go to 1
2      continue
       nat0(n)=m-1
       close(unit=1)
      enddo 


c evaluate the Van der Waals radii
      do n=1,ngroups
       call VdW(60,nat0(n),lab0(n,:),radius(n,:),sym,numb)
      enddo



c read in the Level 3 fragment group numbers
      open(unit=9,file='OUT_FRAGMENTS_LevX',status='old')
      read(9,*)nfbig
      allocate(natbig(nfbig))
      allocate(numbig(nfbig,nsmall))

      do i=1,nfbig
       read(9,*)natbig(i)
       read(9,*)(numbig(i,k),k=1,natbig(i))
      enddo
      close(unit=9)
c read in the Level 1 fragment group numbers
      open(unit=8,file='OUT_FRAGMENTS_Lev1',status='old')
      read(8,*)nffinal
      nffinal1=nffinal
c     write(6,*)' dimension of numat = ',2*nffinal1
      allocate(ksign(2*nffinal1))
      allocate(numat(2*nffinal1))
      allocate(natstore(2*nffinal1,8))
      allocate(nch1(2*nffinal1))
      allocate(nat1(2*nffinal1))
      allocate(c1(2*nffinal1,60,3))
      allocate(lab1(2*nffinal1,60))
      do i=1,nffinal
       read(8,*)numat(i)
       read(8,*)(natstore(i,k),k=1,numat(i))
      enddo
      close(unit=8)
c     stop

c sort natstore (better for cancelling)
      do i=1,nffinal
       call piksrt(numat(i),natstore(i,:))
      enddo
 
c read in the Level 1 signs
      open(unit=7,file='signs.out_Lev1',status='old')
      do k=1,nffinal
       read(7,*)ksign(k)
      enddo
      close(unit=7)

c read the Level1 coordinates
      do n=1,nffinal
       call filelabel(n,ca)
       n1=index(ca,' ') - 1
       ca1='Lev1_COORD'//ca(1:n1)
       open(unit=1,file=ca1,status='old')
       read(1,*)nch1(n)
       m=1
11     continue
       read(1,*,end=21)lab1(n,m),(c1(n,m,k),k=1,3)
       m=m+1
       go to 11
21     continue
       nat1(n)=m-1
       close(unit=1)
      enddo

c the Level 1 fragments may well have repeats, so we cancel these
c save nffinal, as it has been used to dimension arrays
      call cancelL1(2*nffinal1,nffinal,numat,natstore,ksign,nat1,nch1,lab1,c1)

c write out the signs for the L1 frags. Used in the
c induction calculation.
#ifdef __GFORTRAN__
      open(unit=22,file='OUT_L1_FINALSIGNS',status='unknown')
#else
      open(unit=22,file='OUT_L1_FINALSIGNS',status='unknown', 
     . buffered='YES')
#endif
          
      write(22,*)' The number of Lev1 frags after cancelation'
      write(22,*)nffinal
      write(22,*)' The signs'
      do n=1,nffinal
       write(22,*)ksign(n)
      enddo
      close(unit=22)

c nffinal may have been changed
c form exclusion matrix
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

c finished with the Level 3 fragments

c form the distance matrix
      do n1=1,ngroups-1
      do n2=n1+1,ngroups
       call distmin(60,nat0(n1),nat0(n2),c0(n1,:,:),c0(n2,:,:),
     .              radius(n1,:),radius(n2,:),dmin(n1,n2))
       dmin(n2,n1)=dmin(n1,n2)
      enddo
      enddo

      do n1=1,ngroups
       dmin(n1,n1)=0.d0
      enddo


      nffinal0=int(nffinal*nffinal*1.3/2)
      write(6,*)' nffinal = ',nffinal
      write(6,*)' nffinal0 = ',nffinal0
      allocate(matstore(nffinal0,2))
      allocate(nsign(nffinal0))
      do k=1,nffinal0
       nsign(k)=0
       do i=1,2
        matstore(k,i)=0
       enddo
      enddo

      ic=0
      nbextra=0

c start the loops over L1..L1 interactions
      do k1=1,nffinal-1
      do k2=k1+1,nffinal

c     write(6,*)k1,k2,ic

      if(ic.ge.nffinal0-6)then
       write(6,*)' about to overflow'
       write(6,*)k1,k2,ic
       stop
      endif


c find the minimum group-group distance
      dist=1.d6
      do n1=1,numat(k1)
      do n2=1,numat(k2)
       if(dmin(natstore(k1,n1),natstore(k2,n2)).lt.dist)
     .     dist=dmin(natstore(k1,n1),natstore(k2,n2))
      enddo
      enddo


c if the distance is too big, drop this interaction

      jch1=0
      do j1=1,numat(k1)
       jch1=max(jch1,iabs(nch0(natstore(k1,j1))))
      enddo
      jch2=0
      do j2=1,numat(k2)
       jch2=max(jch2,iabs(nch0(natstore(k2,j2))))
      enddo

      nch2=iabs(jch1*jch2)

      if(nch2.eq.0.and.dist.gt.bigtol)go to 111
      if((jch1.eq.0.or.jch2.eq.0).and.dist.gt.bigtol2)go to 111
c if the interaction involves rings go to a special subroutine
      if(numat(k1).gt.2.or.numat(k2).gt.2)then
       ik1k2=ksign(k1)*ksign(k2)
       call intring(ngroups,nom,k1,k2,ik1k2,
     .             ic,matstore,nsign,nffinal0,nffinal1,
     .             numat,natstore,nffinal,nbextra)
       go to 111
      endif
c else proceed for twp groups at most

      do n1=1,2
      do n2=1,2
       mata(n1,n2)=0
      enddo
      enddo

      nsum=0
      do n1=1,numat(k1)
      do n2=1,numat(k2)
       mata(n1,n2)=nom(natstore(k1,n1),natstore(k2,n2))
       nsum=nsum+mata(n1,n2)
      enddo
      enddo

c completely overlapping, drop it
      if(nsum.eq.0)go to 111

      nmax=numat(k1)*numat(k2)



c all allowed
      if(nsum.eq.nmax)then

       ic=ic+1
c record the fragment numbers and the net sign
       matstore(ic,1)=k1
       matstore(ic,2)=k2
       nsign(ic)=ksign(k1)*ksign(k2)
        
       go to 111
      endif


c now partial cases
c use included files for neatness

c (0  1) and (1  0)
      include 'case7.inc'
c (1) and (0)
c (0)     (1) cases
      include 'case8.inc'
c (0  1
c (1  1) and similar cases

c use include files for neatness
      include 'case1.inc'
      include 'case2.inc'
      include 'case3.inc'
      include 'case4.inc'

      if(nsum.eq.2)then

c six cases
        
c (0  0)
c (1  1)
      include 'case9.inc'

c (1  1)
c (0  0)
      include 'case10.inc'

c (1  0)
c (1  0)
      include 'case11.inc'

c (0  1)
c (0  1)
      include 'case12.inc'

c (1  0)
c (0  1)
       include 'case5.inc'
c (0  1)
c (1  0)
       include 'case6.inc'

c end the nsum if
      endif

c final cases where only one element of mata is not zero

c (1  0)
c (0  0) and similar cases
      include 'case13.inc'

111    continue
c end the k1,k2 loops
      enddo
      enddo

      nf=ic

c  write out the final L1 frags
#ifdef __GFORTRAN__
      open(unit=8,file='OUT_FINAL_L1_DATA',status='unknown')
#else
      open(unit=8,file='OUT_FINAL_L1_DATA',status='unknown',
     .           buffered='YES')
#endif
      write(8,*)' Final L1 data'
      write(8,*)' The number of fragments'
      write(8,*)nffinal+nbextra
      write(8,*)' numat and natstore for each'
      do n=1,nffinal+nbextra
       write(8,*)n,numat(n)
       write(8,*)(natstore(n,j),j=1,numat(n))
      enddo
      close(unit=8)

c write out the atoms in the L1 fragments, which will be needed for
c derivative allocations
      nfextra=nffinal+nbextra
      call writeL1atoms(nffinal1,ngroups,nfextra,numat,natstore)

c we now have a complete list of interactions
c and an augmented list of L1 fragments

c we will need to know which groups contain charges
      open(unit=1,file='OUT_CHARGEDGROUPS',status='old')
      read(1,*)
      read(1,*)numcharges
      if(numcharges.gt.0)then
       allocate(numchggrps(numcharges))
       allocate(chchggrps(numcharges))
       read(1,*)
       do k=1,numcharges
        read(1,*)numchggrps(k)
        read(1,*)chchggrps(k)
        read(1,*)
        read(1,*)
       enddo
      endif
      close(unit=1)

c new code 240417 to output groups that provide
c embedded charges for Level=1 (nb) fragments
      if(numcharges.eq.0)go to 2111

      do n=1,nffinal+nbextra
       call filelabel(n,ca1)
       n1=index(ca1,' ')-1
       ca='chnb.'//ca1(1:n1)
#ifdef __GFORTRAN__
       open(unit=21,file=ca,status='unknown')
#else
       open(unit=21,file=ca,status='unknown',buffered='YES')
#endif
       do k=1,numcharges
        do i=1,numat(n)
         if(natstore(n,i).eq.numchggrps(k))go to 2112
        enddo
        write(21,*)numchggrps(k)
2112   enddo
       close(unit=21)
      enddo

2111  continue


      if(numcharges.gt.1)then
c open an output file for long range charge-charge interactions
#ifdef __GFORTRAN__
       open(unit=20,file='OUT_CHARGE_CHARGE',status='unknown')
#else
       open(unit=20,file='OUT_CHARGE_CHARGE',status='unknown',
     .       buffered='YES')
#endif
      endif

c    first we write out the nb.*.0.com files
      call writenbfiles(ngroups,nffinal0,nch0,nat0,lab0,c0,nat1,lab1,c1,
     .        nffinal,nffinal1,nbextra,numcharges,numchggrps,numat,natstore,nch1
     .        ,sym,numb)
c now we print the dalton files
      call printdaltons(ngroups,nat0,nch0,lab0,c0)
c now we write the polar files
      call writepolar(ngroups,nch0,nat0,lab0,c0,sym,numb)

c now the interactions have to be sorted into ab initio and
c electrostatics. The ab initio files need to be written, and
c the details about the remaining interactions ahve to be written
c to a file, so that later the electrostatics can be done.
#ifdef __GFORTRAN__
      open(unit=3,file='OUT_L1L1_data',status='unknown')
      open(unit=33,file='OUT_L1L1_AB_data',status='unknown')
#else
      open(unit=3,file='OUT_L1L1_data',status='unknown',buffered='YES')
      open(unit=33,file='OUT_L1L1_AB_data',status='unknown',
     .            buffered='YES')
#endif

c output the groups involved in ab intio interactions
c so they can later be excluded from dispersion interactions
#ifdef __GFORTRAN__
      open(unit=44,file='OUT_ABGROUPS',status='unknown')
#else
      open(unit=44,file='OUT_ABGROUPS',status='unknown',buffered='YES')
#endif

c write a comment to fort.25 that will be added to OUT_ELECTRONS
      write(25,*)' The number of L1L1 ab calcs is'
      numabcalcs=0
c now do all the interactions

      allocate(n4ab(nfextra*40))
c we will use n4ab to store the valie of n in the loop below for
c ab initio calcs

      do n=1,nf

c identify the L1 fragments
       k1=matstore(n,1)
       k2=matstore(n,2)

c find the minimum frag-frag distance
       dist=1.d6
       do j1=1,numat(k1)
       do j2=1,numat(k2)
        d1=dmin(natstore(k1,j1),natstore(k2,j2))
        if(d1.lt.dist)dist=d1
       enddo
       enddo

c check to see if the L1 fragments contain "charged" groups, that is
c groups that contain formal charges (even if neutral overall)
       if(numcharges.gt.0)then
        match1=0
        do n1=1,numat(k1)
        do k=1,numcharges
c 290119
c        if(natstore(k1,n1).eq.numchggrps(k).
c    .      and.abs(chchggrps(k)).gt.0)match1=1
         if(natstore(k1,n1).eq.numchggrps(k))match1=1
        enddo
        enddo
        match2=0
        do n1=1,numat(k2)
        do k=1,numcharges
c 290119
c         if(natstore(k2,n1).eq.numchggrps(k).
c    .      and.abs(chchggrps(k)).gt.0)match2=1
         if(natstore(k2,n1).eq.numchggrps(k))match2=1
        enddo
        enddo
        if(match1*match2.gt.0)then
c 291019 made xcut = 1.5
        zcut=1.5d0
        if(dtol.gt.zcut)zcut=dtol

c        if(dist.gt.2.d0)then
         if(dist.gt.zcut)then
          write(20,*)(matstore(n,k),k=1,2),nsign(n)
          go to 999
         else
          dist=dtol-1.d0
         endif
        endif
       endif

c 191015, introduce a subroutine that applies only to Level 1 or 2
c The idea is that a L1..L1 interaction cannot be carried out
c ab initio if the caps will clash; rather it will be done
c perturbatively.

       if(Level.eq.2.and.nocapsatall.eq.0)then
        ndim=2*nffinal1
        call Level2exclude(nfextra,ndim,k1,k2,numat,natstore,nflag)
        if(nflag.eq.1)then
         write(45,*)k1,k2,dist
         dist=dtol+1.d0
        endif
       endif

c decide
       if(dist.lt.dtol)then
c its ab initio
        numabcalcs=numabcalcs+1
        n4ab(numabcalcs)=n
       else
c its electrostatics
c all we need to know are the L1 frag numbers and the sign
        write(3,*)(matstore(n,k),k=1,2),nsign(n)
       endif
999   continue
c end the loop over interactions
      enddo
      write(99,*)' before go to 279 ',numabcalcs
      if(numabcalcs.lt.2)go to 279
c some of the ab calcs are repeats, so we have to cancel
      do i=1,numabcalcs-1
       n1=n4ab(i)
       if(nsign(n1).eq.0)go to 4000
       ka1=matstore(n1,1)
       ka2=matstore(n1,2)
       do j=i+1,numabcalcs
        n2=n4ab(j)
        if(nsign(n2).eq.0)go to 4001
        kb1=matstore(n2,1)
        kb2=matstore(n2,2)
        if((ka1.eq.kb1).and.(ka2.eq.kb2))then
         nsign(n1)=nsign(n1)+nsign(n2)
         nsign(n2)=0
        endif
4001   enddo
4000  enddo

      ic=0
      do i=1,numabcalcs
       n1=n4ab(i)
       if(nsign(n1).ne.0)then
        ic=ic+1
        n4ab(ic)=n4ab(i)
       endif
      enddo
      numabcalcs=ic
279   continue

c now write out the files
      if(numabcalcs.eq.0)go to 280
      do j=1,numabcalcs
       n=n4ab(j)
       k1=matstore(n,1)
       k2=matstore(n,2)
        call writeabfiles(ngroups,nffinal0,nf,nch0,nat0,lab0,
     .           c0,nat1,lab1,c1,nffinal,nffinal1,nbextra,numcharges,numchggrps,
     .           numat,natstore,nch1,sym,numb,k1,k2,nsign(n))
        numabgrps=numat(k1)+numat(k2)
        write(44,*)numabgrps
        write(44,*)(natstore(k1,i1),i1=1,numat(k1)),(natstore(k2,i2),i2=1,numat(k2))
        write(33,*)(matstore(n,k),k=1,2),nsign(n)

c new code 240417 to output groups that provide
c embedded charges for Level=1 (nb) fragments
       if(numcharges.eq.0)go to 3111

       call filelabel(k1,ca1)
       call filelabel(k2,ca2)
       n1=index(ca1,' ')-1
       n2=index(ca2,' ')-1
       ca='chab.'//ca1(1:n1)//'.'//ca2(1:n2)
#ifdef __GFORTRAN__
       open(unit=21,file=ca,status='unknown')
#else
       open(unit=21,file=ca,status='unknown',buffered='YES')
#endif
       do k=1,numcharges
        do i=1,numat(k1)
         if(natstore(k1,i).eq.numchggrps(k))go to 3112
        enddo
        do i=1,numat(k2)
         if(natstore(k2,i).eq.numchggrps(k))go to 3112
        enddo
        write(21,*)numchggrps(k)
3112   enddo
       close(unit=21)

3111   continue
      enddo
280   continue
      write(99,*)' after 280 ',numabcalcs
      write(25,*)numabcalcs
      if(numabcalcs.gt.0)then
       write(25,*)' The number of atoms and electrons for each ab calc'
      endif
      close(unit=44)
      close(unit=33)
      close(unit=3)
      close(unit=20)
c we are finished
      end

