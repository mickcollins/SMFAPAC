cmake_minimum_required (VERSION 3.1)
project (SMFAPAC)
enable_language (Fortran)

set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/exe)

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
   # gfortran defaults
   set (CMAKE_Fortran_FLAGS_RELEASE "-cpp -ffixed-line-length-132 -O2")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
   # ifort defaults
   set (CMAKE_Fortran_FLAGS_RELEASE "-132 -O2 -ip -fpe0 -save")
endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")

add_subdirectory(src)
