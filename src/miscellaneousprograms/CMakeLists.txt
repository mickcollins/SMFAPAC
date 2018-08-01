cmake_minimum_required(VERSION 3.1)

add_executable(ABNBderivatives ABNBderivatives.f)
add_executable(adjustcaponly adjustcaponly.f)
add_executable(allocinput allocinput.f)
add_executable(angle angle.f)
add_executable(bond bond.f)
add_executable(combineallderivatives combineallderivatives.f)
add_executable(convert2dalton convert2dalton.f)
add_executable(convertcoords4games convertcoords4games.f)
add_executable(convertgamhessian convertgamhessian.f)
add_executable(dihedral dihedral.f)
add_executable(distinctlabels distinctlabels.f)
add_executable(findHandCbonds_VdW findHandCbonds_VdW.f)
add_executable(gamstone2cart gamstone2cart.f)
add_executable(makeIN_CONSTRAINTS makeIN_CONSTRAINTS.f)
add_executable(newHatom newHatom.f)
add_executable(orderlist orderlist.f)
add_executable(setprocs setprocs.f)
add_executable(Stonepunch2cart_Matt Stonepunch2cart_Matt.f)
add_executable(updatehessian updatehessian.f)
add_executable(PAA PAA.f)

install (TARGETS 
   ABNBderivatives
   adjustcaponly
   allocinput
   angle
   bond
   combineallderivatives
   convert2dalton
   convertcoords4games
   convertgamhessian
   dihedral
   distinctlabels
   findHandCbonds_VdW
   gamstone2cart
   makeIN_CONSTRAINTS
   newHatom
   orderlist
   setprocs
   Stonepunch2cart_Matt
   updatehessian
   PAA
 DESTINATION ${CMAKE_INSTALL_PREFIX})
