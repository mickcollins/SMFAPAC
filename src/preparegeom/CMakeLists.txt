cmake_minimum_required(VERSION 3.1)

set(EXE Preparegeom)

set(SRC
   AtomicData.f90
   adjustsingle.f
   anum2elneg.f
   assignbyelect.f
   bondlimits.f
   elements.f
#   functions.f  # seems to have duplicate function definitions
   indexi.f
   indexx.f
   specgrp.f
)

add_executable( ${EXE} ${SRC})
install (TARGETS ${EXE} DESTINATION ${CMAKE_INSTALL_PREFIX})
