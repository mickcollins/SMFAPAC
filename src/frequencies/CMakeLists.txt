cmake_minimum_required(VERSION 3.1)

set(EXE Freq ${SRC})

set(SRC
   freq.f
   project.f
)

add_executable( ${EXE} ${SRC})
install (TARGETS ${EXE} DESTINATION ${CMAKE_INSTALL_PREFIX})

find_package(BLAS)
if (BLAS_FOUND)
   target_link_libraries(${EXE} ${BLAS_LIBRARIES})
endif()
