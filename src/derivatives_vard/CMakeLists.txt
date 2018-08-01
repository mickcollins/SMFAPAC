cmake_minimum_required(VERSION 3.1)

set(EXE Derivs_smfa ${SRC})

set(SRC
   anum2alab.f
   derivcomb.f
   derivheader.f90
   massset.f
   readdata.f
   readfrag.f
   readsigns.f
)

add_executable( ${EXE} ${SRC})
install (TARGETS ${EXE} DESTINATION ${CMAKE_INSTALL_PREFIX})
