cmake_minimum_required(VERSION 3.1)

set(EXE gdma)

set(SRC
   atom_grids.f90
   dma.f90
   gdma.f90
   input.F90
   timing.F90
)

add_executable( ${EXE} ${SRC})
install (TARGETS ${EXE} DESTINATION ${CMAKE_INSTALL_PREFIX})
