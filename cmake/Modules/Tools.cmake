if(BUILD_TOOLS)
  add_executable(binary2txt ${LAMMPS_TOOLS_DIR}/binary2txt.cpp)
  target_compile_definitions(binary2txt PRIVATE -DLAMMPS_${LAMMPS_SIZES})
  install(TARGETS binary2txt DESTINATION ${CMAKE_INSTALL_BINDIR})

  include(CheckGeneratorSupport)
  if(CMAKE_GENERATOR_SUPPORT_FORTRAN)
    include(CheckLanguage)
    check_language(Fortran)
    if(CMAKE_Fortran_COMPILER)
      enable_language(Fortran)
      add_executable(chain.x ${LAMMPS_TOOLS_DIR}/chain.f90)
      target_link_libraries(chain.x PRIVATE ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES})
      add_executable(micelle2d.x ${LAMMPS_TOOLS_DIR}/micelle2d.f90)
      target_link_libraries(micelle2d.x PRIVATE ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES})
      install(TARGETS chain.x micelle2d.x DESTINATION ${CMAKE_INSTALL_BINDIR})
    else()
      message(WARNING "No suitable Fortran compiler found, skipping build of 'chain.x' and 'micelle2d.x'")
    endif()
  else()
    message(WARNING "CMake build doesn't support Fortran, skipping build of 'chain.x' and 'micelle2d.x'")
  endif()

  enable_language(C)
  get_filename_component(MSI2LMP_SOURCE_DIR ${LAMMPS_TOOLS_DIR}/msi2lmp/src ABSOLUTE)
  file(GLOB MSI2LMP_SOURCES ${MSI2LMP_SOURCE_DIR}/[^.]*.c)
  add_executable(msi2lmp ${MSI2LMP_SOURCES})
  if(STANDARD_MATH_LIB)
    target_link_libraries(msi2lmp PRIVATE ${STANDARD_MATH_LIB})
  endif()
  install(TARGETS msi2lmp DESTINATION ${CMAKE_INSTALL_BINDIR})
  install(FILES ${LAMMPS_DOC_DIR}/msi2lmp.1 DESTINATION ${CMAKE_INSTALL_MANDIR}/man1)
endif()

if(BUILD_LAMMPS_SHELL)
  find_package(PkgConfig REQUIRED)
  pkg_check_modules(READLINE IMPORTED_TARGET REQUIRED readline)
  if(NOT LAMMPS_EXCEPTIONS)
    message(WARNING "The LAMMPS shell needs LAMMPS_EXCEPTIONS enabled for full functionality")
  endif()

  # include resource compiler to embed icons into the executable on Windows
  if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
    enable_language(RC)
    set(ICON_RC_FILE ${LAMMPS_TOOLS_DIR}/lammps-shell/lmpicons.rc)
  endif()

  add_executable(lammps-shell ${LAMMPS_TOOLS_DIR}/lammps-shell/lammps-shell.cpp ${ICON_RC_FILE})
  target_include_directories(lammps-shell PRIVATE ${LAMMPS_TOOLS_DIR}/lammps-shell)

  # workaround for broken readline pkg-config file on FreeBSD
  if(CMAKE_SYSTEM_NAME STREQUAL "FreeBSD")
    target_include_directories(lammps-shell PRIVATE /usr/local/include)
  endif()
  target_link_libraries(lammps-shell PRIVATE lammps PkgConfig::READLINE)
  install(TARGETS lammps-shell EXPORT LAMMPS_Targets DESTINATION ${CMAKE_INSTALL_BINDIR})
  install(DIRECTORY ${LAMMPS_TOOLS_DIR}/lammps-shell/icons DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/)
  install(FILES ${LAMMPS_TOOLS_DIR}/lammps-shell/lammps-shell.desktop DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/applications/)
endif()


