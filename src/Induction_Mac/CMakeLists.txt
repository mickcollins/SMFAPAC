cmake_minimum_required(VERSION 3.1)

set( EXE Induct_Mac)

set(SRC
   Tfactor2.f
   atomicnumber.f
   filelabel.f
   fld.f
   induction.f
   readfamilies.f
   readgroupdata.f
   readmoments.f
)

add_executable( ${EXE} ${SRC})
install (TARGETS ${EXE} DESTINATION ${CMAKE_INSTALL_PREFIX})
