cmake_minimum_required(VERSION 3.1)

set(EXE Disp)

set(SRC
   Disp.f
)

add_executable( ${EXE} ${SRC})
install (TARGETS ${EXE} DESTINATION ${CMAKE_INSTALL_PREFIX})
