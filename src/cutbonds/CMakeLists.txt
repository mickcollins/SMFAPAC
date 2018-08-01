cmake_minimum_required(VERSION 3.1)

set(EXE cutsmfa)

set(SRC
   bigfract.f
   breakgroup.f
   cancel.f
   checkbonds.f
   compress.f
   eight.f
   expand.f
   extrasubs.f
   families.f
   filelabel.f
   finalcancel.f
   findfragelim.f
   five.f
   four.f
   fractdata.f90
   fragsteps.f
   makefrag.f
   piksrt.f
   readINFRACT.f
   readin_big.f
   setupfrags1.f
   seven.f
   shrink.f
   six.f
   three.f
   writecom.f
   writefrags.f
   writerawdata.f
)

add_executable( ${EXE} ${SRC})
install (TARGETS ${EXE} DESTINATION ${CMAKE_INSTALL_PREFIX})
