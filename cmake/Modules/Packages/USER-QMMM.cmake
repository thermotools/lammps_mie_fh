enable_language(C)

if(NOT BUILD_SHARED_LIBS)
  message(WARNING "It is recommended to use BUILD_SHARED_LIBS=yes with USER-QMMM")
endif()
add_library(qmmm STATIC ${LAMMPS_LIB_SOURCE_DIR}/qmmm/libqmmm.c)
if(NOT BUILD_SHARED_LIBS)
  install(TARGETS qmmm EXPORT LAMMPS_Targets LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()
set_target_properties(qmmm PROPERTIES OUTPUT_NAME lammps_qmmm${LAMMPS_LIB_SUFFIX})
target_link_libraries(lammps PRIVATE qmmm)
target_include_directories(qmmm PUBLIC ${LAMMPS_LIB_SOURCE_DIR}/qmmm)
