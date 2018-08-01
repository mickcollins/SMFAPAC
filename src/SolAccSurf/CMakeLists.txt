cmake_minimum_required(VERSION 3.1)

add_executable( SolAccSurface SolAccSurface.f)
add_executable( makecube makecube.f )

install (TARGETS SolAccSurface makecube DESTINATION ${CMAKE_INSTALL_PREFIX})
