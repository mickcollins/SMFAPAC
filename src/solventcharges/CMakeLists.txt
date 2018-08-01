cmake_minimum_required(VERSION 3.1)

set(EXE SOLCH)

set(SRC
   checkbonds.f
   fractdata.f90
   main.f
   readin_big.f
   solcharges.f
)

add_executable( ${EXE} ${SRC})
install (TARGETS ${EXE} DESTINATION ${CMAKE_INSTALL_PREFIX})
