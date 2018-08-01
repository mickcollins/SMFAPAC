      subroutine massset
      use derivheader

      integer, allocatable  :: numa(:)

       real*8 atomic_masses(110)

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


      allocate(numa(natom))

c get the atomic number and coordinates for each real atom
      do m=1,nfrag
      do n=1,nat0(m)
      numa(nat(m,n))=numstore(m,n)
      do k=1,3
      coord(nat(m,n),k)=c(m,n,k)
      enddo
      enddo
      enddo

c replace coordinates with the higher precision values
c contained in name.xyz
      open(unit=50,file='name.xyz',status='old')
      read(50,*)
      read(50,*)
      do n=1,natom
       read(50,*)lab(n),(coord(n,k),k=1,3)
      enddo
      close(unit=50)

      bohr=1.d0/1.8897259886d0

      do k=1,3
      do n=1,natom
       coord(n,k)=coord(n,k)/bohr
      enddo
      enddo

c find the mass and element symbol for each real atom
      do n=1,natom
       na=numa(n)
       ams=atomic_masses(na)
       amas(n)=ams
       call anum2alab(na,lab(n))
      enddo


      return
      end
