if(PKG_COMPRESS)
  find_package(ZLIB REQUIRED)
  include_directories(${ZLIB_INCLUDE_DIRS})
  list(APPEND LAMMPS_LINK_LIBS ${ZLIB_LIBRARIES})
endif()
