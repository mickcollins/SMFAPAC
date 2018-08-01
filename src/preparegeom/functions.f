      real*8 function anum2cov(na)


c in Angstrom, from Mick's atomic.dat, other values converted from Chris
c Mick's atomic.dat misses 85-87 incl and from 95 onward

      integer :: na
      real*8 :: covalent_radii(110)

      data covalent_radii/
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


      anum2cov=covalent_radii(na)

      return

      end function anum2cov

      integer function alab2anum(alab)

c given the atom label this returns the atomic number

      character*2 :: alab,atomic_labels(110)

      integer n

      data atomic_labels/'H','He','Li','Be','B','C','N','O','F','Ne',
     .'Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr',
     .'Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb',
     . 'Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',
     . 'Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu',
     . 'Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os',
     . 'Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac',
     . 'Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No',
     . 'Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds'/

      do n=1,110
       if(atomic_labels(n).eq.alab)then
        alab2anum=n
        return
       endif
      enddo
      alab2anum=0

      return

      end function alab2anum

      real*8 function anum2elneg(na)

c the electronegativities of the elements
c in atomic number order

      integer :: na
      real*8 :: electroneg(110)

      data electroneg/
     . 2.300,4.160,0.912,1.576,2.051,2.544,3.066,3.610,4.193,4.789,
     . 0.869,1.293,1.613,1.916,2.253,2.589,2.869,3.242,0.734,1.034,
     . 1.19,1.38,1.53,1.65,1.75,1.80,1.84,1.88,1.85,1.59,1.756,1.994,
     . 2.211,2.434,2.685,2.966,0.706,0.963,1.12,1.32,1.41,1.47,1.51,
     . 1.54,1.56,1.59,1.87,1.52,1.656,1.824,1.984,2.158,2.359,2.582,
     . 0.659,0.881,
     . 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
     . 1.09,1.16,1.34,1.47,1.60,1.65,1.68,1.72,1.92,1.76,
     . 1.789,1.854,2.01,2.19,2.39,2.60,0.67,0.89,
     . 0.,0.,0.,0.,0.,0.,
     . 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0./

      anum2elneg=electroneg(na)

      return

      end function anum2elneg

