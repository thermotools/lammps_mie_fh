set(WIN_PACKAGES ASPHERE BODY CLASS2 COLLOID COMPRESS CORESHELL DIPOLE GPU
                 GRANULAR KSPACE LATTE MANYBODY MC MISC ML-IAP MOLECULE OPT
                 PERI POEMS QEQ REPLICA RIGID SHOCK ML-SNAP SPIN SRD VORONOI
                 ATC AWPMD BOCS BROWNIAN CG-DNA CG-SDK
                 COLVARS DIFFRACTION DPD-REACT DRUDE EFF FEP
                 ML-HDNNP INTEL MANIFOLD MDI MEAM DPD-MESO
                 MESONT USER-MISC MGPT MOFFF MOLFILE OPENMP
                 PHONON PTM QTB REACTION REAXFF
                 DPD-SMOOTH MACHDYN SMTBQ SPH TALLY UEF
                 YAFF DIELECTRIC)

foreach(PKG ${WIN_PACKAGES})
  set(PKG_${PKG} ON CACHE BOOL "" FORCE)
endforeach()

# these two packages require a full MPI implementation
if(BUILD_MPI)
  set(PKG_MPIIO ON CACHE BOOL "" FORCE)
  set(PKG_LATBOLTZ ON CACHE BOOL "" FORCE)
endif()

set(DOWNLOAD_VORO ON CACHE BOOL "" FORCE)
set(DOWNLOAD_EIGEN3 ON CACHE BOOL "" FORCE)
set(LAMMPS_MEMALIGN "0" CACHE STRING "" FORCE)
set(CMAKE_TUNE_FLAGS "-Wno-missing-include-dirs" CACHE STRING "" FORCE)
set(CMAKE_EXE_LINKER_FLAGS "-Wl,--enable-stdcall-fixup,--as-needed,-lssp" CACHE STRING "" FORCE)
set(CMAKE_SHARED_LINKER_FLAGS "-Wl,--enable-stdcall-fixup,--as-needed,-lssp" CACHE STRING "" FORCE)
set(BUILD_TOOLS ON CACHE BOOL "" FORCE)
set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_BINARY_DIR}/lammps-installer")
