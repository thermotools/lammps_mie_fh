# Preset that turns on all existing packages. Using the combination
# of this preset followed by the nolib.cmake preset should configure
# a LAMMPS binary, with as many packages included, that can be compiled
# with just a working C++ compiler and an MPI library.

set(ALL_PACKAGES
  ADIOS
  ASPHERE
  ATC
  AWPMD
  BOCS
  BODY
  BROWNIAN
  CG-DNA
  CG-SDK
  CLASS2
  COLLOID
  COLVARS
  COMPRESS
  CORESHELL
  DIELECTRIC
  DIFFRACTION
  DIPOLE
  DPD-BASIC
  DPD-MESO
  DPD-REACT
  DPD-SMOOTH
  DRUDE
  EFF
  FEP
  GPU
  GRANULAR
  H5MD
  INTEL
  INTERLAYER
  KIM
  KOKKOS
  KSPACE
  LATBOLTZ
  LATTE
  MACHDYN
  MANIFOLD
  MANYBODY
  MC
  MDI
  MEAM
  MESONT
  MESSAGE
  MGPT
  MISC
  ML-HDNNP
  ML-IAP
  ML-PACE
  ML-QUIP
  ML-RANN
  ML-SNAP
  MOFFF
  MOLECULE
  MOLFILE
  MPIIO
  MSCG
  NETCDF
  OPENMP
  OPT
  PERI
  PHONON
  PLUGIN
  PLUMED
  POEMS
  PTM
  PYTHON
  QEQ
  QMMM
  QTB
  REACTION
  REAXFF
  REPLICA
  RIGID
  SCAFACOS
  SHOCK
  SMTBQ
  SPH
  SPIN
  SRD
  TALLY
  UEF
  USER-MISC
  VORONOI
  VTK
  YAFF)

foreach(PKG ${ALL_PACKAGES})
  set(PKG_${PKG} ON CACHE BOOL "" FORCE)
endforeach()
