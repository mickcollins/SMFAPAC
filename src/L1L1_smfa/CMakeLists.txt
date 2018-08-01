cmake_minimum_required(VERSION 3.1)

set(EXE L1L1_SMFA)

set(SRC
   AtNum.f
   L1L1_Mac.f
   Level2exclude.f
   VdW.f
   anum2vdw.f
   cancelL1.f
   cancelL1L1.f
   dalton_files_Mac.f90
   distmin.f
   filelabel.f
   intring.f
   match_mg.f
   molecule_adt.f90
   piksrt.f
   printdaltons.f
   writeL1atoms.f
   writeabfiles.f
   writenbfiles.f
   writepolar.f
)

add_executable( ${EXE} ${SRC})
install (TARGETS ${EXE} DESTINATION ${CMAKE_INSTALL_PREFIX})
