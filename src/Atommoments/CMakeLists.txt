cmake_minimum_required(VERSION 3.1)

set(EXE AllMoments)

set(SRC
   filelabel.f
   momentcomb.f
   readdata.f
   readfrag.f
   readsigns.f
   momentheader.f90
)

add_executable( ${EXE} ${SRC})
install (TARGETS ${EXE} DESTINATION ${CMAKE_INSTALL_PREFIX})
