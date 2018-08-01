      subroutine solcharges

      use fractdata

      real*8, allocatable :: fammass(:)

      integer, allocatable  :: natid(:),natch(:)

      real*8 atomic_masses(110)

      character*2 sollab(40)

      integer alab2anum

       data atomic_masses/
     .    1.0079,   4.0026,   6.9410,   9.0122,  10.8110,  12.0107,
     .   14.0067,  15.9994,  18.9984,  20.1797,  22.9898,  24.3050,
     .   26.9815,  28.0855,  30.9738,  32.0650,  35.4530,  39.9480,
     .   39.0983,  40.0780,  44.9559,  47.8670,  50.9415,  51.9961,
     .   54.9380,  55.8450,  58.6934,  58.9332,  63.5460,  65.3800,
     .   69.7230,  72.6400,  74.9216,  78.9600,  79.9040,  83.7980,
     .   85.4678,  87.6200,  88.9059,  91.2240,  92.9064,  95.9600,
     .   98.0000, 101.0700, 102.9055, 106.4200, 107.8682, 112.4110,
     .  114.8180, 118.7100, 121.7600, 127.6000, 126.9045, 131.2930,
     .  132.9055, 137.3270, 138.9055, 140.1160, 140.9077, 144.2420,
     .  145.0000, 150.3600, 151.9640, 157.2500, 158.9254, 162.5000,
     .  164.9303, 167.2590, 168.9342, 173.0540, 174.9668, 178.4900,
     .  180.9479, 183.8400, 186.2070, 190.2300, 192.2170, 195.0840,
     .  196.9666, 200.5900, 204.3833, 207.2000, 208.9804, 210.0000,
     .  210.0000, 220.0000, 223.0000, 226.0000, 227.0000, 231.0359,
     .  232.0381, 237.0000, 238.0289, 243.0000, 244.0000, 247.0000,
     .  247.0000, 251.0000, 252.0000, 257.0000, 258.0000, 259.0000,
     .  262.0000, 261.0000, 262.0000, 266.0000, 264.0000, 277.0000,
     .  268.0000, 271.0000/


      allocate(fammass(natom))

c calculate the mass of each family

      do n=1,natom
       fammass(n)=0.d0
       do i=1,nfam(n)
        num=numa(ifam(n,i))
        fammass(n)=fammass(n)+atomic_masses(num)
       enddo
      enddo

c we expect that fammass identities the chemical composition exactly

c read in the chemical composition of the solvent molecule
c that we want to be included in "IN_CHARGES"

      open(unit=1,file="IN_SOLCHARGES",status='unknown')
      read(1,*,end=2000)
      read(1,*)nsol
      read(1,*)
      solmass=0.d0
      do i=1,nsol
       read(1,*)sollab(i)
       nblank=index(sollab(i)," ")
       if(nblank.eq.2)then
        num=alab2anum(sollab(i)(1:1))
       else
        num=alab2anum(sollab(i))
       endif
       solmass=solmass+atomic_masses(num)

      enddo
      close(unit=1)
100   format(a2)



c now each family that has the same mass as solmass is to be "charged"

      newcharges=0
      do n=1,natom
       if(abs(solmass-fammass(n)).lt.1.d-2)then
        newcharges=newcharges+1
       endif
      enddo

c get the original versin of IN_CHARGES
      numold=0
      open(unit=1,file='IN_CHARGES',status='unknown')
      read(1,*,end=1)
      read(1,*)numold
1      ntotal=numold+newcharges
       allocate(natid(ntotal))
       allocate(natch(ntotal))
      if(numold.gt.0)then
       read(1,*)
       do i=1,numold
        read(1,*)natid(i),natch(i)
       enddo
      endif
      close(unit=1)

c add the new charges
      n1=numold
      do n=1,natom
       if(abs(solmass-fammass(n)).lt.1.d-2)then
        n1=n1+1
        natid(n1)=ifam(n,1)
        natch(n1)=0
       endif
      enddo

c output the new version
      open(unit=1,file='IN_CHARGES',status='unknown')
      write(1,*)" The number of charged groups is"
      write(1,*)ntotal
      write(1,*) " The atom numbers and charges are"
      do n=1,ntotal
       write(1,*)natid(n),natch(n)
      enddo
      close(unit=1)

2000  return
      end
