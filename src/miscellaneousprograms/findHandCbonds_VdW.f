      program findHandC_VdW

c this program finds hydrogen bonds, which can be included in the
c bonding description, with a multiplicity of -1
c In the fragmentation process, when bonds with multiplicity of -1
c are broken, no capping H atom is introduced.

c The program also defines a bond with multiplicity of -1
c between groups (families) that are both charged and close together


      implicit double precision(a-h,o-z)
      parameter(maxfam=100)

      real*8, allocatable  :: c(:,:),radius(:),amas(:)
      real*8, allocatable  :: dist(:,:),distg(:,:),avdw(:,:)
      real*8, allocatable  :: disthbonds(:)
      real*8, allocatable :: VdWrad(:)

      integer, allocatable :: numa(:),ndb(:),mdb(:),link(:,:)
      integer, allocatable :: ichg(:)
      integer, allocatable :: nHbonds(:,:),nprint(:)
      integer, allocatable :: natstore(:,:),nfam(:),mfam(:)
      integer, allocatable :: ifam(:,:),ngpch(:)
      integer, allocatable :: linkg(:,:)
      integer, allocatable :: nval(:),nhbondedto(:)

      character*2, allocatable :: lab(:)


      character*2 sym(110)
      dimension rad(110), am(110)
      dimension numb(110)


      character*5 Htrue, Ctrue

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


      Bohr=1.d0/1.88972598860

c see what has to be done

      Htrue='true'
      Ctrue='true'
      factor=1.0d0
      open(unit=3,file='IN_MINUSBONDS',status='unknown')
      read(3,*,end=201)
      read(3,200)Htrue
      read(3,*)
      read(3,200)Ctrue
      read(3,*)
      read(3,*)factor
c factor multiplies the sum of the VdW radii
200   format(a5)
      close(unit=3)
      if(Htrue(1:4).eq.'true'.or.Htrue(1:5).eq.'false')go to 201
      if(Ctrue(1:4).eq.'true'.or.Ctrue(1:5).eq.'false')go to 201
      write(6,*)' only true or false allows in IN_MINUSBONDS'
      stop
201   continue

      hdist=2.4d0

c get the number of atoms from name.mol

      open(unit=1,file='name.mol',status='old')
      read(1,*)
      read(1,*)
      read(1,*)
      read(1,*)natom,nbonds

      allocate(lab(natom))
      allocate(c(natom,3))
      allocate(radius(natom))
      allocate(amas(natom))
      allocate(numa(natom))
      allocate(ndb(4*natom))
      allocate(mdb(4*natom))
      allocate(link(natom,natom))
      allocate(dist(natom,natom))
      allocate(ichg(natom))
      allocate(nHbonds(2*natom,2))
      allocate(nprint(2*natom))
      allocate(disthbonds(2*natom))
      allocate(natstore(2*natom,natom))

      allocate(VdWrad(natom))
      allocate(nval(natom))
      allocate(nhbondedto(natom))

      do i=1,natom
      read(1,8001)(c(i,k),k=1,3),lab(i),ichg(i)
      enddo
8001  format(3f10.4,1x,a2,i4)

      do n=1,natom
      do m=1,natom
       link(n,m)=0
      enddo
      enddo

c define the properties of each atom
      do n=1,natom
      match=0
      do m=1,110
        if(lab(n).eq.sym(m))then
        numa(n)=m
        amas(n)=am(m)
        radius(n)=rad(m)
        match=1
        endif
      enddo
      if(match.eq.0)then
       write(6,*)' The atom number ',n,' is not recognised'
       write(6,*)' by atomic number or element symbol'
       stop
      endif

      enddo

c assign VdW radii
      do n=1,natom
       VdWrad(n)=anum2vdw(numa(n))*bohr
      enddo
c assign the number of valence electrons 
c valence is used in hydrogen bonding
      do n=1,natom
       nval(n)=0
       if(numa(n).le.2)nval(n)=numa(n)
       if(numa(n).gt.2.and.numa(n).le.10)nval(n)=numa(n)-2
       if(numa(n).gt.10.and.numa(n).le.18)nval(n)=numa(n)-10
       if(numa(n).gt.18.and.numa(n).le.20)nval(n)=numa(n)-18
       if(numa(n).gt.30.and.numa(n).le.36)nval(n)=numa(n)-28
       if(numa(n).gt.36.and.numa(n).le.38)nval(n)=numa(n)-36
       if(numa(n).gt.48.and.numa(n).le.54)nval(n)=numa(n)-46
       if(numa(n).gt.54.and.numa(n).le.56)nval(n)=numa(n)-54
       if(numa(n).gt.80.and.numa(n).le.86)nval(n)=numa(n)-78
      enddo

c zero nhbondedto array
      do n=1,natom
       nhbondedto(n)=0
      enddo


c read in the families
      open(unit=20,file='families.out',status='old')

      read(20,*)
      read(20,*)ngroup
c allocate memory
      allocate(nfam(ngroup))
      allocate(mfam(natom))
      allocate(ifam(ngroup,maxfam))
      allocate(ngpch(ngroup))
      allocate(distg(ngroup,ngroup))
      allocate(avdw(ngroup,ngroup))
      allocate(linkg(ngroup,ngroup))
      do n1=1,ngroup
      do n2=1,ngroup
       linkg(n1,n2)=0
      enddo
      enddo

      read(20,*)
      do n=1,ngroup
       read(20,*)nfam(n)
c nfam is the number of atom in each group
       read(20,*)(ifam(n,i),i=1,nfam(n))
c ifam is the identity of the atoms in each group
      enddo
      close(unit=20)

c assign each family to each atom
       do n=1,ngroup
       do i=1,nfam(n)
        mfam(ifam(n,i))=n
       enddo
       enddo

c now read the bonds from name.mol
c linkg = 1 if two groups are connected by a real bond
      do i=1,nbonds
       read(1,*)n1,m1
       linkg(mfam(n1),mfam(m1))=1
      enddo
 
      close(unit=1)

c  calculate the bondlengths
      do n=1,natom-1
      do m=n+1,natom
      sum=0.d0
      do k=1,3
      sum=sum+(c(m,k)-c(n,k))**2
      enddo
      sum=sqrt(sum)
      dist(n,m)=sum
      dist(m,n)=sum
c     write(6,*)n,m,sum
c fix the tolerance for this pair
c and mark bonds
c link=1 if two atoms are connected by a real bond
      tol=radius(n)+radius(m)+0.4d0
      if(sum.lt.tol)then
c     write(69,*)n,m,sum
      link(n,m)=1
      link(m,n)=1
      endif

      enddo
      enddo

c make sure a H atom is only linked to the closeast heavy atom
      do n=1,natom
       if(numa(n).eq.1)then
        ic=0
        iclose=0
        distmax=1.d6
        do m=1,natom
         if(link(n,m).eq.1)then
          ic=ic+1
          if(dist(n,m).lt.distmax)then
           iclose=m
           distmax=dist(n,m)
          endif
         endif
        enddo
        if(ic.ge.1)nhbondedto(n)=iclose
       endif
      enddo

          
c the cos of the Hbond angle must be less than a given value
c labelled climit. Taken to be the value for 135 degrees.
      climit=-0.7071d0

c now look for H bonds
c new code trying to have same effect as that in prepareliquid
      ic=0
      if(Htrue(1:4).eq.'true')then
       do n=1,natom
        if(nhbondedto(n).ne.0)then
         m=nhbondedto(n)
         if(nval(m).ge.5.and.nval(m).le.7)then
          do k=1,natom
           if(nval(k).ge.5.and.nval(k).le.7)then
            tol=VdWrad(n)+VdWrad(k)-0.32d0
            if(dist(n,k).le.tol)then
c should check the angle of the hydrogen bond
c and only allow "near" linear
c first calc the angle
             rnm=0.d0
             rnk=0.d0
             dot=0.d0
             do j=1,3
              rnm=rnm+(c(m,j)-c(n,j))**2
              rnk=rnk+(c(k,j)-c(n,j))**2
              dot=dot+(c(m,j)-c(n,j))*(c(k,j)-c(n,j))
             enddo
             cang=dot/sqrt(rnm*rnk)
             if(cang.lt.climit)then
              ic=ic+1
              nhbonds(ic,1)=m
              nhbonds(ic,2)=k
              disthbonds(ic)=dist(n,k)
             endif
            endif
           endif
          enddo
         endif
        endif
       enddo
      endif

c      ic=0
c      if(Htrue(1:4).eq.'true')then
c
c      do n=1,natom
c       if(lab(n).ne.'H ')go to 21
c      do m=1,natom
c       if(lab(m).eq.'C ')go to 22
c       if(link(n,m).eq.0)go to 22
c       do k=1,natom
c        if(link(n,k).eq.1)go to 23
c        if(m.eq.k)go to 23
c        if(lab(k).eq.'C ')go to 23
c        if(lab(k).eq.'H ')go to 23
c        if(dist(n,k).le.hdist)then
c        ic=ic+1
c        nhbonds(ic,1)=m
c        nhbonds(ic,2)=k
c        disthbonds(ic)=dist(n,k)
cc       write(69,*)n,m,k
c        endif
c23      enddo
c22      enddo
c21      enddo

c      endif


c now look for "charge" bonds
c skip if "false"
      if(Ctrue.eq.'false')go to 2000

c check for no charges
      ncharges=0
      do i=1,natom
      if(ichg(i).ne.0)ncharges=ncharges+1
      enddo
      if(ncharges.eq.0)go to 2000


c assign charges to families (if any atom is charged)
      do n=1,ngroup
       ngpch(n)=0
       do j=1,nfam(n)
        if(ichg(ifam(n,j)).ne.0)ngpch(n)=1
       enddo
      enddo

c calculate the minimum group-group distances relative to
c the vdw radii
      do n=1,ngroup-1
      do m=n+1,ngroup
      ratmin=1.d6
      do i=1,nfam(n)
      do j=1,nfam(m)
       svdW=anum2vdw(numa(ifam(n,i)))+anum2vdw(numa(ifam(m,j)))
      sum=0.d0
      do k=1,3
      sum=sum+(c(ifam(n,i),k)-c(ifam(m,j),k))**2
      enddo
       ratio=sqrt(sum)/svdW
      if(ratio.lt.ratmin)then
         ratmin=ratio
         dmin=sum
         vdW=svdW
       endif
      enddo
      enddo
      distg(n,m)=sqrt(dmin)
      distg(m,n)=distg(n,m)
      avdw(n,m)=vdW
      avdw(m,n)=vdW
      enddo
      enddo

c assign "C bonds" between the first atoms in groups
      do n=1,ngroup-1
      do m=n+1,ngroup
       if(ngpch(n)*ngpch(m).eq.0)go to 3001
c convert Van der Waals radii to Angstrom from Bohr
       dn2m=avdw(n,m)*factor*0.529177d0
       if(distg(n,m).gt.dn2m)go to 3001
c we create a new H-bond
       ic=ic+1
       nhbonds(ic,1)=ifam(n,1)
       nhbonds(ic,2)=ifam(m,1)
       disthbonds(ic)=distg(n,m)
3001  enddo
      enddo

2000  continue


c check for and remove repeats

      if(ic.eq.1)then
       nic=1
       nprint(1)=1
      endif

      if(ic.gt.1)then

      do n=1,ic
       nprint(n)=1
      enddo

c do not allow H or C bonds within a group
      do n=1,ic-1
      do m=n+1,ic
       if(mfam(nhbonds(n,1)).eq.mfam(nhbonds(m,1)).and.
     &    mfam(nhbonds(n,2)).eq.mfam(nhbonds(m,2)))nprint(m)=0
       if(mfam(nhbonds(n,1)).eq.mfam(nhbonds(m,2)).and.
     &    mfam(nhbonds(n,2)).eq.mfam(nhbonds(m,1)))nprint(m)=0

      enddo
      enddo

c checking
c     write(69,*)'first nprint'
c     do i=1,ic
c      write(69,*)nhbonds(i,1),nhbonds(i,2),nprint(i)
c     enddo



c do not allow H or C bonds between groups that are already
c connected by real bonds
      do n=1,ic
       if(linkg(mfam(nhbonds(n,1)),mfam(nhbonds(n,2))).eq.1)nprint(n)=0
      enddo

      endif

c checking
c     write(69,*)'second nprint'
c     do i=1,ic
c      write(69,*)nhbonds(i,1),nhbonds(i,2),nprint(i)
c     enddo

       nbd=ic

      nic=0
      do n=1,ic
      if(nprint(n).eq.1)nic=nic+1
      enddo



      if(ic.gt.0)then
      open(unit=7,file='IN_HBONDS',status='unknown')
      write(7,*)' The number of Hydrogen bonds is'
      write(7,*)nic
      write(7,*)' The heavy atoms connected by these are '
      do i=1,ic
      if(nprint(i).eq.1)then
      write(7,*)nhbonds(i,1),nhbonds(i,2),-1
      endif
      enddo
      close(unit=7)
      endif

      end

      real*8 function anum2vdw(na)

      integer na
      real*8 vanderwaals_radii(110)

      data vanderwaals_radii/ 
     &2.2676727, 2.6456187, 3.4393037, 3.7794547, 3.7794547, 3.2125367,
     &2.9290777, 2.8723857, 2.7778997, 2.9101797, 4.2896807, 3.2692277,
     &3.7794547, 3.9684267, 3.4015087, 3.4015087, 3.3070227, 3.5526877,
     &5.1967497, 3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.7794547,
     &3.7794547, 3.7794547, 3.7794547, 3.0802557, 2.6456187, 2.6267207,
     &3.5337897, 3.7794547, 3.4959957, 3.5904817, 3.4959957, 3.8172487,
     &3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.7794547,
     &3.7794547, 3.7794547, 3.7794547, 3.0802557, 3.2503307, 2.9857687,
     &3.6471737, 4.1007077, 3.7794547, 3.8928377, 3.7416597, 4.0818107,
     &3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.7794547,
     &3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.7794547,
     &3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.7794547,
     &3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.2503307,
     &3.1369477, 2.9290777, 3.7038657, 3.8172487, 3.7794547, 3.7794547,
     &3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.7794547,
     &3.7794547, 3.5148927, 3.7794547, 3.7794547, 3.7794547, 3.7794547,
     &3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.7794547,
     &3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.7794547,
     &3.7794547, 3.7794547/

       anum2vdw=vanderwaals_radii(na)

       return

       end function anum2vdw

