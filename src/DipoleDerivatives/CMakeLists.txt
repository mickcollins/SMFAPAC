cmake_minimum_required(VERSION 3.1)

set(EXE DipoleDerivs)

set(SRC
   derivcomb.f
   derivheader.f90
   readdata.f
   readfrag.f
   readsigns.f
)

add_executable( ${EXE} ${SRC})
install (TARGETS ${EXE} DESTINATION ${CMAKE_INSTALL_PREFIX})
