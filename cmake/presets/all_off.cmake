# preset that turns on all existing packages off. can be used to reset
# an existing package selection without losing any other settings

set(ALL_PACKAGES ASPHERE BODY CLASS2 COLLOID COMPRESS CORESHELL DIPOLE GPU
        GRANULAR KIM KOKKOS KSPACE LATTE MANYBODY MC MISC MESSAGE MOLECULE
        MPIIO MSCG OPT PERI POEMS PYTHON QEQ REPLICA RIGID SHOCK SNAP SPIN
        SRD VORONOI
        USER-ADIOS USER-ATC USER-AWPMD USER-BOCS USER-CGDNA USER-CGSDK
        USER-COLVARS USER-DIFFRACTION USER-DPD USER-DRUDE USER-EFF USER-FEP
        USER-H5MD USER-INTEL USER-LB USER-MANIFOLD USER-MEAMC USER-MESODPD
        USER-MGPT USER-MISC USER-MOFFF USER-MOLFILE USER-NETCDF USER-OMP
        USER-PHONON USER-PLUMED USER-PTM USER-QMMM USER-QTB USER-QUIP
        USER-REACTION USER-REAXC USER-SCAFACOS USER-SDPD USER-SMD USER-SMTBQ
        USER-SPH USER-TALLY USER-UEF USER-VTK USER-YAFF)

foreach(PKG ${ALL_PACKAGES})
  set(PKG_${PKG} OFF CACHE BOOL "" FORCE)
endforeach()
