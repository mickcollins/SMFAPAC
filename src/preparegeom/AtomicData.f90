module AtomicData
 implicit none

 contains

integer function aNum(i)
! a function which returns the atomic number corresponding to its entry
! index.
! Note that after Polonium (element 84) atomic number does not match
! index number 
 implicit none
 integer::i
 integer, dimension(91)::AtomicNumbers

! initialise the data
 data AtomicNumbers / 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, &
 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, &
 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, &
 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, &
 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 88, 89, 90, 91, 92, 93, 94/

 aNum = AtomicNumbers(i)

 return

end function aNum

character(len=2) function aLab(i)
 implicit none
! a function which returns the atomic label corresponding to its entry
! index.
! Note that after Polonium (element 84) atomic number does not match
! index number 
 integer::i
 character(len=2), dimension(91)::AtomicLabel

!initialise the data
 data AtomicLabel / "H ", "He", "Li", "Be", "B ", "C ", "N ", "O ", "F", &
 "Ne", "Na", "Mg", "Al", "Si", "P ", "S ", "Cl", "Ar", "K ", "Ca", "Sc", &
 "Ti", "V ", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", &
 "Se", "Br", "Kr", "Rb", "Sr", "Y ", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", &
 "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I ", "Xe", "Cs", "Ba", "La", &
 "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", &
 "Yb", "Lu", "Hf", "Ta", "W ", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", &
 "Pb", "Bi", "Po", "Ra", "Ac", "Th", "Pa", "U ", "Np", "Pu"/

 aLab=AtomicLabel(i)

 return

end function aLab

real function aMass(i)
 implicit none
! a function which returns the atomic mass corresponding to its entry
! index.
! Note that after Polonium (element 84) atomic number does not match
! index number 
 integer::i
 real, dimension(91)::AtomicMass

! initialise the data
 data AtomicMass/ 1.00794, 4.002602, 6.941, 9.012182, 10.811, 12.0107, &
 14.0067, 15.9994, 18.9984032, 20.1797, 22.98977, 24.305, 26.981538, &
 28.0855, 30.973761, 32.065, 35.453, 39.948, 39.0983, 40.078, 44.95591,&
 47.867, 50.9415, 51.9961, 54.938049, 55.845, 58.9332, 58.6934, 63.546,&
 65.409, 69.723, 72.64, 74.9216, 78.96, 79.904, 83.798, 85.4678, 87.62,&
 88.90585, 91.224, 92.90638, 95.94, 98., 101.07, 102.9055, 106.42, &
107.8682, &
 112.411, 114.818, 118.71, 121.76, 127.6, 126.90447, 131.293, 132.90545, &
 137.327, 138.9055, 140.116, 140.90765, 144.24, 145., 150.36, 151.964, & 
 157.25, 158.92534, 162.5, 164.93032, 167.259, 168.93421, 173.04,      &
174.967, &
 178.49, 180.9479, 183.84, 186.207, 190.23, 192.217, 195.078, 196.96655,&
 200.59, 204.3833, 207.2, 208.98038, 209., 226., 227., 232.0381, &
231.03588, &
 238.02891, 237., 244./

 aMass=AtomicMass(i)

 return

end function aMass

real function aRad(i)
! a function which returns the atomic radii corresponding to its entry
! index.
! Note that after Polonium (element 84) atomic number does not match
! index number 
 integer::i
 real, dimension(91)::AtomicRadii

! initialise the data
 data AtomicRadii / 0.23, 1.22, 0.68, 0.35, 0.83, 0.68, 0.68, 0.68,    &
0.64, 1.6, &
 0.97, 1.1, 1.35, 1.2, 1.05, 1.02, 0.99, 1.92, 1.33, 0.99, 1.44, 1.47, &
1.33, &
 0.67, 1.35, 1.34, 1.33, 1.5, 1.52, 1.45, 1.22, 1.17, 1.21, 1.22, 1.21, &
1.98, &
 1.47, 1.12, 1.78, 1.56, 1.48, 1.47, 1.35, 1.4, 1.45, 1.5, 1.59, 1.69, &
1.63, &
 1.46, 1.46, 1.47, 1.4, 2.18, 1.67, 1.34, 1.87, 1.83, 1.82, 1.81, 1.8, &
1.8, &
 1.99, 1.79, 1.76, 1.75, 1.74, 1.73, 1.72, 1.94, 1.72, 1.57, 1.43, 1.37, &
1.35, &
 1.37, 1.32, 1.5, 1.5, 1.7, 1.55, 1.54, 1.54, 1.68, 1.9, 1.88, 1.79, &
1.61, &
 1.58, 1.55, 1.53/

 aRad=AtomicRadii(i)

 return

end function aRad

integer function aLabToANum(label)
! a function that takes an atomic label as input and returns the
! corresponding atomic number
 character(len=2), intent(in)::label
 
 integer::i 
 logical::match
 
!sort through the atomic labels to find a match
 match=.false.
 do i=1, 91
  if (aLab(i) == label) then
   match=.true.
   aLabToANum=aNum(i)
  end if
 end do
 if (match .eqv. .false.) then
  print *, 'In function aLabToANum: unrecognised atomic label'
  stop 1
 end if

 return

end function aLabToANum

logical function CheckAtomValid(label)
 implicit none
! A function that checks whether the atomic label supplied as an
! argument appear in the file atomic.dat and hence whether it will be
! recognised by other algorithms written
 character(len=2), intent(in)::label

 integer::i
 logical::valid

 valid=.false.
 do i=1, 91
  if(label == aLab(i)) then
   valid=.true.
   exit
  end if
 end do

 CheckAtomValid=valid

 return

end function CheckAtomValid

real function ANumToVDW(atomicNumber)
 implicit none
! a small program the converts atomic number to a vanderwaal radii
 real, parameter::bohr=0.52917721
 
 integer, intent(in)::atomicNumber

 real, dimension(110)::vanderwaals_radii

 data vanderwaals_radii/ &
  2.2676727, 2.6456187, 3.4393037, 3.7794547, 3.7794547, 3.2125367, &
  2.9290777, 2.8723857, 2.7778997, 2.9101797, 4.2896807, 3.2692277, &
  3.7794547, 3.9684267, 3.4015087, 3.4015087, 3.3070227, 3.5526877, &
  5.1967497, 3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.7794547, &
  3.7794547, 3.7794547, 3.7794547, 3.0802557, 2.6456187, 2.6267207, &
  3.5337897, 3.7794547, 3.4959957, 3.5904817, 3.4959957, 3.8172487, &
  3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.7794547, &
  3.7794547, 3.7794547, 3.7794547, 3.0802557, 3.2503307, 2.9857687, &
  3.6471737, 4.1007077, 3.7794547, 3.8928377, 3.7416597, 4.0818107, &
  3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.7794547, &
  3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.7794547, &
  3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.7794547, &
  3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.2503307, &
  3.1369477, 2.9290777, 3.7038657, 3.8172487, 3.7794547, 3.7794547, &
  3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.7794547, &
  3.7794547, 3.5148927, 3.7794547, 3.7794547, 3.7794547, 3.7794547, &
  3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.7794547, &
  3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.7794547, 3.7794547, &
  3.7794547, 3.7794547/

 ANumToVDW=vanderwaals_radii(atomicNumber)*bohr

 return

end function ANumToVDW


end module AtomicData


