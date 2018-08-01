cmake_minimum_required(VERSION 3.1)

set( EXE Elect_Mac )

set(SRC
   Tfactorfast.f
   ct0c_d.f
   ct1d.f
   ct1d_d.f
   ct2q.f
   ct2q_d.f
   ct3o.f
   ct3o_d.f
   ct4h.f
   ct4h_d.f
   dt2d.f
   dt2d_d.f
   dt3q.f
   dt3q_d.f
   dt4o.f
   dt4o_d.f
   elect.f
   filelabel.f
   nbelect_Mac.f
   qt4q.f
   qt4q_d.f
   symm.f
   zero.f
)

add_executable( ${EXE} ${SRC})
install (TARGETS ${EXE} DESTINATION ${CMAKE_INSTALL_PREFIX})
