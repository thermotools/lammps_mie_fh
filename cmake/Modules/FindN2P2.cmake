include(FindPackageHandleStandardArgs)

# Check if N2P2_DIR is set manually.
if (DEFINED ENV{N2P2_DIR})
  set(N2P2_DIR "${N2P2_DIR}")
# If not, try if directory "lib/nnp/n2p2" exists.
else()
  get_filename_component(_fullpath "${LAMMPS_LIB_SOURCE_DIR}/nnp/n2p2" REALPATH)
  if (EXISTS ${_fullpath})
    set(N2P2_DIR "${_fullpath}")
  endif()
endif()

# Set path to include directory.
find_path(N2P2_INCLUDE_DIR InterfaceLammps.h PATHS "${N2P2_DIR}/include")
# Two libraries need to be linked: libnnp and libnnpif.
find_library(N2P2_LIBNNP NAMES nnp PATHS "${N2P2_DIR}/lib")
find_library(N2P2_LIBNNPIF NAMES nnpif PATHS "${N2P2_DIR}/lib")
# Users could compile n2p2 with extra flags which are then also required for
# pair_nnp.cpp compilation. To forward them to the LAMMPS build process n2p2
# writes a file with cmake commands, e.g.
#
# target_compile_definitions(lammps PRIVATE -DNNP_NO_SF_GROUPS)
#
# to "lib/lammps-extra.cmake" which is then included by USER-NNP.cmake.
find_file(N2P2_CMAKE_EXTRA NAMES lammps-extra.cmake PATHS "${N2P2_DIR}/lib")

find_package_handle_standard_args(N2P2 DEFAULT_MSG
  N2P2_DIR
  N2P2_INCLUDE_DIR
  N2P2_LIBNNP
  N2P2_LIBNNPIF
  N2P2_CMAKE_EXTRA)

if(N2P2_FOUND)
  set(N2P2_INCLUDE_DIRS ${N2P2_INCLUDE_DIR})
  set(N2P2_LIBRARIES ${N2P2_LIBNNPIF} ${N2P2_LIBNNP})
  set(N2P2_CMAKE_EXTRAS ${N2P2_CMAKE_EXTRA})

  mark_as_advanced(
    N2P2_DIR
    N2P2_INCLUDE_DIR
    N2P2_LIBNNP
    N2P2_LIBNNPIF
    N2P2_CMAKE_EXTRA
  )
endif()
