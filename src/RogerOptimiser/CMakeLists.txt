cmake_minimum_required(VERSION 3.1)

set(EXE ROGEROPT)

set(SRC
   constraints.f
   dangle.f
   dbond.f
   ddihed.f
   hesaugTS.f
   optimiser.f
   tsvector.f
)


add_executable( ${EXE} ${SRC})
install (TARGETS ${EXE} DESTINATION ${CMAKE_INSTALL_PREFIX})

find_package(BLAS)
if (BLAS_FOUND)
   target_link_libraries(${EXE} ${BLAS_LIBRARIES})
endif()
