      program convertcoords4games
      implicit double precision(a-h,o-z)

      character*2 lab
      real*8 c(3)
      real an
      integer alab2anum

c     open(unit=1,file='gamesscoords',status='unknown')

1     continue
      read(5,*,end=2)lab,(c(k),k=1,3)
      n=alab2anum(lab)
      an=float(n)
c     write(1,100)lab,an,(c(k),k=1,3)
      write(6,100)lab,an,(c(k),k=1,3)
100   format(a2,f6.1,3f13.6)
      go to 1
2     continue
c     close(unit=1)
      end


      integer function alab2anum(alab)

c given the atom label this returns the atomic number

      character*2 :: alab,atomic_labels(110)

      data atomic_labels/'H','He','Li','Be','B','C','N','O','F','Ne',
     .'Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr',
     . 'Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb',
     . 'Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',
     . 'Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu',
     . 'Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os',
     . 'Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac',
     .   'Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No'
     . ,'Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds'/

      integer n

      do n=1,110
        if(atomic_labels(n).eq.alab)then
      alab2anum=n
      return
        endif
      enddo
      alab2anum=0

      return
      end
