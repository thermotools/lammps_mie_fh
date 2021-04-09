set(MOLFILE_INCLUDE_DIR "${LAMMPS_LIB_SOURCE_DIR}/molfile" CACHE STRING "Path to VMD molfile plugin headers")
set(MOLFILE_INCLUDE_DIRS "${MOLFILE_INCLUDE_DIR}")
add_library(molfile INTERFACE)
target_include_directories(molfile INTERFACE ${MOLFILE_INCLUDE_DIRS})
target_link_libraries(lammps PRIVATE molfile)
