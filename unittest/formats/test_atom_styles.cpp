/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "atom.h"
#include "input.h"
#include "lammps.h"
#include "utils.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cstdio>
#include <cstring>
#include <mpi.h>

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

namespace LAMMPS_NS {
using ::testing::Eq;

class AtomStyleTest : public ::testing::Test {
protected:
    LAMMPS *lmp;

    void SetUp() override
    {
        const char *args[] = {"SimpleCommandsTest", "-log", "none", "-echo", "screen", "-nocite"};
        char **argv        = (char **)args;
        int argc           = sizeof(args) / sizeof(char *);
        if (!verbose) ::testing::internal::CaptureStdout();
        lmp = new LAMMPS(argc, argv, MPI_COMM_WORLD);
        if (!verbose) ::testing::internal::GetCapturedStdout();
        ASSERT_NE(lmp, nullptr);
        if (!verbose) ::testing::internal::CaptureStdout();
        lmp->input->one("units real");
        lmp->input->one("dimension 3");
        lmp->input->one("pair_style zero 4.0");
        lmp->input->one("region box block -4 4 -4 4 -4 4");
        if (!verbose) ::testing::internal::GetCapturedStdout();
    }

    void TearDown() override
    {
        if (!verbose) ::testing::internal::CaptureStdout();
        delete lmp;
        if (!verbose) ::testing::internal::GetCapturedStdout();
        remove("test_atom_styles.data");
        remove("test_atom_styles.restart");
    }
};

TEST_F(AtomStyleTest, atomic)
{
    ASSERT_THAT(std::string(lmp->atom->atom_style), Eq("atomic"));
    ASSERT_NE(lmp->atom->avec, nullptr);
    ASSERT_EQ(lmp->atom->natoms, 0);
    ASSERT_EQ(lmp->atom->nlocal, 0);
    ASSERT_EQ(lmp->atom->nghost, 0);
    ASSERT_EQ(lmp->atom->nmax, 1);
    ASSERT_EQ(lmp->atom->tag_enable, 1);
    ASSERT_EQ(lmp->atom->molecular, 0);
    ASSERT_EQ(lmp->atom->nellipsoids, 0);
    ASSERT_EQ(lmp->atom->nlines, 0);
    ASSERT_EQ(lmp->atom->ntris, 0);
    ASSERT_EQ(lmp->atom->nbodies, 0);
    ASSERT_EQ(lmp->atom->nbonds, 0);
    ASSERT_EQ(lmp->atom->nangles, 0);
    ASSERT_EQ(lmp->atom->ndihedrals, 0);
    ASSERT_EQ(lmp->atom->nimpropers, 0);
    ASSERT_EQ(lmp->atom->ntypes, 0);
    ASSERT_EQ(lmp->atom->nbondtypes, 0);
    ASSERT_EQ(lmp->atom->nangletypes, 0);
    ASSERT_EQ(lmp->atom->ndihedraltypes, 0);
    ASSERT_EQ(lmp->atom->nimpropertypes, 0);
    ASSERT_EQ(lmp->atom->bond_per_atom, 0);
    ASSERT_EQ(lmp->atom->angle_per_atom, 0);
    ASSERT_EQ(lmp->atom->dihedral_per_atom, 0);
    ASSERT_EQ(lmp->atom->improper_per_atom, 0);
    ASSERT_EQ(lmp->atom->extra_bond_per_atom, 0);
    ASSERT_EQ(lmp->atom->extra_angle_per_atom, 0);
    ASSERT_EQ(lmp->atom->extra_dihedral_per_atom, 0);
    ASSERT_EQ(lmp->atom->extra_improper_per_atom, 0);

    ASSERT_EQ(lmp->atom->sphere_flag, 0);
    ASSERT_EQ(lmp->atom->ellipsoid_flag, 0);
    ASSERT_EQ(lmp->atom->line_flag, 0);
    ASSERT_EQ(lmp->atom->tri_flag, 0);
    ASSERT_EQ(lmp->atom->body_flag, 0);
    ASSERT_EQ(lmp->atom->peri_flag, 0);
    ASSERT_EQ(lmp->atom->electron_flag, 0);
    ASSERT_EQ(lmp->atom->wavepacket_flag, 0);
    ASSERT_EQ(lmp->atom->sph_flag, 0);
    ASSERT_EQ(lmp->atom->molecule_flag, 0);
    ASSERT_EQ(lmp->atom->molindex_flag, 0);
    ASSERT_EQ(lmp->atom->molatom_flag, 0);
    ASSERT_EQ(lmp->atom->q_flag, 0);
    ASSERT_EQ(lmp->atom->mu_flag, 0);
    ASSERT_EQ(lmp->atom->rmass_flag, 0);
    ASSERT_EQ(lmp->atom->radius_flag, 0);
    ASSERT_EQ(lmp->atom->omega_flag, 0);
    ASSERT_EQ(lmp->atom->torque_flag, 0);
    ASSERT_EQ(lmp->atom->angmom_flag, 0);
    ASSERT_EQ(lmp->atom->vfrac_flag, 0);
    ASSERT_EQ(lmp->atom->spin_flag, 0);
    ASSERT_EQ(lmp->atom->eradius_flag, 0);
    ASSERT_EQ(lmp->atom->ervel_flag, 0);
    ASSERT_EQ(lmp->atom->erforce_flag, 0);
    ASSERT_EQ(lmp->atom->cs_flag, 0);
    ASSERT_EQ(lmp->atom->csforce_flag, 0);
    ASSERT_EQ(lmp->atom->vforce_flag, 0);
    ASSERT_EQ(lmp->atom->ervelforce_flag, 0);
    ASSERT_EQ(lmp->atom->etag_flag, 0);
    ASSERT_EQ(lmp->atom->rho_flag, 0);
    ASSERT_EQ(lmp->atom->esph_flag, 0);
    ASSERT_EQ(lmp->atom->cv_flag, 0);
    ASSERT_EQ(lmp->atom->vest_flag, 0);
    ASSERT_EQ(lmp->atom->dpd_flag, 0);
    ASSERT_EQ(lmp->atom->edpd_flag, 0);
    ASSERT_EQ(lmp->atom->tdpd_flag, 0);
    ASSERT_EQ(lmp->atom->mesont_flag, 0);
    ASSERT_EQ(lmp->atom->sp_flag, 0);
    ASSERT_EQ(lmp->atom->x0_flag, 0);
    ASSERT_EQ(lmp->atom->smd_flag, 0);
    ASSERT_EQ(lmp->atom->damage_flag, 0);
    ASSERT_EQ(lmp->atom->contact_radius_flag, 0);
    ASSERT_EQ(lmp->atom->smd_data_9_flag, 0);
    ASSERT_EQ(lmp->atom->smd_stress_flag, 0);
    ASSERT_EQ(lmp->atom->eff_plastic_strain_flag, 0);
    ASSERT_EQ(lmp->atom->eff_plastic_strain_rate_flag, 0);
    ASSERT_EQ(lmp->atom->pdscale, 1.0);

    ASSERT_NE(lmp->atom->tag, nullptr);
    ASSERT_NE(lmp->atom->type, nullptr);
    ASSERT_NE(lmp->atom->mask, nullptr);
    ASSERT_NE(lmp->atom->image, nullptr);
    ASSERT_NE(lmp->atom->x, nullptr);
    ASSERT_NE(lmp->atom->v, nullptr);
    ASSERT_NE(lmp->atom->f, nullptr);
    ASSERT_EQ(lmp->atom->q, nullptr);
    ASSERT_EQ(lmp->atom->mu, nullptr);
    ASSERT_EQ(lmp->atom->omega, nullptr);
    ASSERT_EQ(lmp->atom->angmom, nullptr);
    ASSERT_EQ(lmp->atom->torque, nullptr);
    ASSERT_EQ(lmp->atom->radius, nullptr);
    ASSERT_EQ(lmp->atom->rmass, nullptr);
    ASSERT_EQ(lmp->atom->ellipsoid, nullptr);
    ASSERT_EQ(lmp->atom->line, nullptr);
    ASSERT_EQ(lmp->atom->tri, nullptr);
    ASSERT_EQ(lmp->atom->body, nullptr);
    ASSERT_EQ(lmp->atom->molecule, nullptr);
    ASSERT_EQ(lmp->atom->molindex, nullptr);
    ASSERT_EQ(lmp->atom->molatom, nullptr);
    ASSERT_EQ(lmp->atom->num_bond, nullptr);
    ASSERT_EQ(lmp->atom->bond_type, nullptr);
    ASSERT_EQ(lmp->atom->bond_atom, nullptr);
    ASSERT_EQ(lmp->atom->num_angle, nullptr);
    ASSERT_EQ(lmp->atom->angle_type, nullptr);
    ASSERT_EQ(lmp->atom->angle_atom1, nullptr);
    ASSERT_EQ(lmp->atom->angle_atom2, nullptr);
    ASSERT_EQ(lmp->atom->angle_atom3, nullptr);
    ASSERT_EQ(lmp->atom->num_dihedral, nullptr);
    ASSERT_EQ(lmp->atom->dihedral_type, nullptr);
    ASSERT_EQ(lmp->atom->dihedral_atom1, nullptr);
    ASSERT_EQ(lmp->atom->dihedral_atom2, nullptr);
    ASSERT_EQ(lmp->atom->dihedral_atom3, nullptr);
    ASSERT_EQ(lmp->atom->dihedral_atom4, nullptr);
    ASSERT_EQ(lmp->atom->num_improper, nullptr);
    ASSERT_EQ(lmp->atom->improper_type, nullptr);
    ASSERT_EQ(lmp->atom->improper_atom1, nullptr);
    ASSERT_EQ(lmp->atom->improper_atom2, nullptr);
    ASSERT_EQ(lmp->atom->improper_atom3, nullptr);
    ASSERT_EQ(lmp->atom->improper_atom4, nullptr);
    ASSERT_EQ(lmp->atom->maxspecial, 1);
    ASSERT_EQ(lmp->atom->nspecial, nullptr);
    ASSERT_EQ(lmp->atom->special, nullptr);
    ASSERT_EQ(lmp->atom->vfrac, nullptr);
    ASSERT_EQ(lmp->atom->s0, nullptr);
    ASSERT_EQ(lmp->atom->x0, nullptr);
    ASSERT_EQ(lmp->atom->sp, nullptr);
    ASSERT_EQ(lmp->atom->fm, nullptr);
    ASSERT_EQ(lmp->atom->fm_long, nullptr);
    ASSERT_EQ(lmp->atom->spin, nullptr);
    ASSERT_EQ(lmp->atom->eradius, nullptr);
    ASSERT_EQ(lmp->atom->ervel, nullptr);
    ASSERT_EQ(lmp->atom->erforce, nullptr);
    ASSERT_EQ(lmp->atom->ervelforce, nullptr);
    ASSERT_EQ(lmp->atom->cs, nullptr);
    ASSERT_EQ(lmp->atom->csforce, nullptr);
    ASSERT_EQ(lmp->atom->vforce, nullptr);
    ASSERT_EQ(lmp->atom->etag, nullptr);
    ASSERT_EQ(lmp->atom->uCond, nullptr);
    ASSERT_EQ(lmp->atom->uMech, nullptr);
    ASSERT_EQ(lmp->atom->uChem, nullptr);
    ASSERT_EQ(lmp->atom->uCG, nullptr);
    ASSERT_EQ(lmp->atom->uCGnew, nullptr);
    ASSERT_EQ(lmp->atom->duChem, nullptr);
    ASSERT_EQ(lmp->atom->dpdTheta, nullptr);
    ASSERT_EQ(lmp->atom->cc, nullptr);
    ASSERT_EQ(lmp->atom->cc_flux, nullptr);
    ASSERT_EQ(lmp->atom->edpd_temp, nullptr);
    ASSERT_EQ(lmp->atom->edpd_flux, nullptr);
    ASSERT_EQ(lmp->atom->edpd_cv, nullptr);
    ASSERT_EQ(lmp->atom->length, nullptr);
    ASSERT_EQ(lmp->atom->buckling, nullptr);
    ASSERT_EQ(lmp->atom->bond_nt, nullptr);
    ASSERT_EQ(lmp->atom->contact_radius, nullptr);
    ASSERT_EQ(lmp->atom->smd_data_9, nullptr);
    ASSERT_EQ(lmp->atom->smd_stress, nullptr);
    ASSERT_EQ(lmp->atom->eff_plastic_strain, nullptr);
    ASSERT_EQ(lmp->atom->eff_plastic_strain_rate, nullptr);
    ASSERT_EQ(lmp->atom->damage, nullptr);
    ASSERT_EQ(lmp->atom->rho, nullptr);
    ASSERT_EQ(lmp->atom->drho, nullptr);
    ASSERT_EQ(lmp->atom->esph, nullptr);
    ASSERT_EQ(lmp->atom->desph, nullptr);
    ASSERT_EQ(lmp->atom->cv, nullptr);
    ASSERT_EQ(lmp->atom->vest, nullptr);
    ASSERT_EQ(lmp->atom->nmolecule, 0);
    ASSERT_EQ(lmp->atom->molecules, nullptr);
    ASSERT_EQ(lmp->atom->nivector, 0);
    ASSERT_EQ(lmp->atom->ndvector, 0);
    ASSERT_EQ(lmp->atom->iname, nullptr);
    ASSERT_EQ(lmp->atom->dname, nullptr);
    ASSERT_EQ(lmp->atom->mass, nullptr);
    ASSERT_EQ(lmp->atom->mass_setflag, nullptr);
    ASSERT_EQ(lmp->atom->nextra_grow, 0);
    ASSERT_EQ(lmp->atom->nextra_restart, 0);
    ASSERT_EQ(lmp->atom->nextra_border, 0);
    ASSERT_EQ(lmp->atom->nextra_grow_max, 0);
    ASSERT_EQ(lmp->atom->nextra_restart_max, 0);
    ASSERT_EQ(lmp->atom->nextra_border_max, 0);
    ASSERT_EQ(lmp->atom->nextra_store, 0);
    ASSERT_EQ(lmp->atom->extra_grow, nullptr);
    ASSERT_EQ(lmp->atom->extra_restart, nullptr);
    ASSERT_EQ(lmp->atom->extra_border, nullptr);
    ASSERT_EQ(lmp->atom->extra, nullptr);
    ASSERT_EQ(lmp->atom->sametag, nullptr);
    ASSERT_EQ(lmp->atom->map_style, 0);
    ASSERT_EQ(lmp->atom->map_user, 0);
    ASSERT_EQ(lmp->atom->map_tag_max, -1);

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("atom_style charge");
    lmp->input->one("atom_style atomic");
    if (!verbose) ::testing::internal::GetCapturedStdout();

    ASSERT_THAT(std::string(lmp->atom->atom_style), Eq("atomic"));
    ASSERT_NE(lmp->atom->avec, nullptr);
    ASSERT_EQ(lmp->atom->natoms, 0);
    ASSERT_EQ(lmp->atom->nlocal, 0);
    ASSERT_EQ(lmp->atom->nghost, 0);
    ASSERT_EQ(lmp->atom->nmax, 1);
    ASSERT_EQ(lmp->atom->tag_enable, 1);
    ASSERT_EQ(lmp->atom->molecular, 0);
    ASSERT_EQ(lmp->atom->ntypes, 0);

    ASSERT_EQ(lmp->atom->molecule_flag, 0);
    ASSERT_EQ(lmp->atom->molindex_flag, 0);
    ASSERT_EQ(lmp->atom->molatom_flag, 0);

    ASSERT_EQ(lmp->atom->q_flag, 0);
    ASSERT_EQ(lmp->atom->q, nullptr);

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("atom_modify map hash");
    lmp->input->one("create_box 2 box");
    lmp->input->one("create_atoms 1 single -2.0  2.0  0.1");
    lmp->input->one("create_atoms 1 single -2.0 -2.0 -0.1");
    lmp->input->one("create_atoms 2 single  2.0  2.0 -0.1");
    lmp->input->one("create_atoms 2 single  2.0 -2.0  0.1");
    lmp->input->one("mass 1 4.0");
    lmp->input->one("mass 2 2.4");
    lmp->input->one("pair_coeff * *");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_THAT(std::string(lmp->atom->atom_style), Eq("atomic"));
    ASSERT_NE(lmp->atom->avec, nullptr);
    ASSERT_EQ(lmp->atom->natoms, 4);
    ASSERT_EQ(lmp->atom->nlocal, 4);
    ASSERT_EQ(lmp->atom->nghost, 0);
    ASSERT_NE(lmp->atom->nmax, -1);
    ASSERT_EQ(lmp->atom->tag_enable, 1);
    ASSERT_EQ(lmp->atom->molecular, 0);
    ASSERT_EQ(lmp->atom->ntypes, 2);

    ASSERT_NE(lmp->atom->mass, nullptr);
    ASSERT_NE(lmp->atom->mass_setflag, nullptr);
    ASSERT_NE(lmp->atom->sametag, nullptr);
    ASSERT_EQ(lmp->atom->map_style, 2);
    ASSERT_EQ(lmp->atom->map_user, 2);
    ASSERT_EQ(lmp->atom->map_tag_max, 4);
    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("pair_coeff * *");
    lmp->input->one("write_data test_atom_styles.data nocoeff");
    lmp->input->one("clear");
    lmp->input->one("atom_style atomic");
    lmp->input->one("pair_style zero 4.0");
    lmp->input->one("atom_modify map array");
    lmp->input->one("units real");
    lmp->input->one("read_data test_atom_styles.data");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_THAT(std::string(lmp->atom->atom_style), Eq("atomic"));
    ASSERT_NE(lmp->atom->avec, nullptr);
    ASSERT_EQ(lmp->atom->natoms, 4);
    ASSERT_EQ(lmp->atom->nlocal, 4);
    ASSERT_EQ(lmp->atom->nghost, 0);
    ASSERT_NE(lmp->atom->nmax, -1);
    ASSERT_EQ(lmp->atom->tag_enable, 1);
    ASSERT_EQ(lmp->atom->molecular, 0);
    ASSERT_EQ(lmp->atom->ntypes, 2);

    double **x = lmp->atom->x;
    double **v = lmp->atom->v;
    ASSERT_DOUBLE_EQ(x[0][0], -2.0);
    ASSERT_DOUBLE_EQ(x[0][1], 2.0);
    ASSERT_DOUBLE_EQ(x[0][2], 0.1);
    ASSERT_DOUBLE_EQ(x[1][0], -2.0);
    ASSERT_DOUBLE_EQ(x[1][1], -2.0);
    ASSERT_DOUBLE_EQ(x[1][2], -0.1);
    ASSERT_DOUBLE_EQ(x[2][0], 2.0);
    ASSERT_DOUBLE_EQ(x[2][1], 2.0);
    ASSERT_DOUBLE_EQ(x[2][2], -0.1);
    ASSERT_DOUBLE_EQ(x[3][0], 2.0);
    ASSERT_DOUBLE_EQ(x[3][1], -2.0);
    ASSERT_DOUBLE_EQ(x[3][2], 0.1);
    ASSERT_DOUBLE_EQ(v[0][0], 0.0);
    ASSERT_DOUBLE_EQ(v[0][1], 0.0);
    ASSERT_DOUBLE_EQ(v[0][2], 0.0);
    ASSERT_DOUBLE_EQ(v[1][0], 0.0);
    ASSERT_DOUBLE_EQ(v[1][1], 0.0);
    ASSERT_DOUBLE_EQ(v[1][2], 0.0);
    ASSERT_DOUBLE_EQ(v[2][0], 0.0);
    ASSERT_DOUBLE_EQ(v[2][1], 0.0);
    ASSERT_DOUBLE_EQ(v[2][2], 0.0);
    ASSERT_DOUBLE_EQ(v[3][0], 0.0);
    ASSERT_DOUBLE_EQ(v[3][1], 0.0);
    ASSERT_DOUBLE_EQ(v[3][2], 0.0);

    ASSERT_DOUBLE_EQ(lmp->atom->mass[1], 4.0);
    ASSERT_DOUBLE_EQ(lmp->atom->mass[2], 2.4);
    ASSERT_EQ(lmp->atom->mass_setflag[1], 1);
    ASSERT_EQ(lmp->atom->mass_setflag[2], 1);
    ASSERT_EQ(lmp->atom->map_style, 1);
    ASSERT_EQ(lmp->atom->map_user, 1);
    ASSERT_EQ(lmp->atom->map_tag_max, 4);
    ASSERT_EQ(lmp->atom->tag_consecutive(), 1);

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("pair_coeff * *");
    lmp->input->one("group two id 2:4:2");
    lmp->input->one("delete_atoms group two compress no");
    lmp->input->one("write_restart test_atom_styles.restart");
    lmp->input->one("clear");
    lmp->input->one("read_restart test_atom_styles.restart");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_THAT(std::string(lmp->atom->atom_style), Eq("atomic"));
    ASSERT_NE(lmp->atom->avec, nullptr);
    ASSERT_EQ(lmp->atom->natoms, 2);
    ASSERT_EQ(lmp->atom->nlocal, 2);
    ASSERT_EQ(lmp->atom->nghost, 0);
    ASSERT_NE(lmp->atom->nmax, -1);
    ASSERT_EQ(lmp->atom->tag_enable, 1);
    ASSERT_EQ(lmp->atom->molecular, 0);
    ASSERT_EQ(lmp->atom->ntypes, 2);
    ASSERT_EQ(lmp->atom->tag_consecutive(), 0);

    x = lmp->atom->x;
    v = lmp->atom->v;
    ASSERT_DOUBLE_EQ(x[0][0], -2.0);
    ASSERT_DOUBLE_EQ(x[0][1], 2.0);
    ASSERT_DOUBLE_EQ(x[0][2], 0.1);
    ASSERT_DOUBLE_EQ(x[1][0], 2.0);
    ASSERT_DOUBLE_EQ(x[1][1], 2.0);
    ASSERT_DOUBLE_EQ(x[1][2], -0.1);
    ASSERT_DOUBLE_EQ(v[0][0], 0.0);
    ASSERT_DOUBLE_EQ(v[0][1], 0.0);
    ASSERT_DOUBLE_EQ(v[0][2], 0.0);
    ASSERT_DOUBLE_EQ(v[1][0], 0.0);
    ASSERT_DOUBLE_EQ(v[1][1], 0.0);
    ASSERT_DOUBLE_EQ(v[1][2], 0.0);

    ASSERT_DOUBLE_EQ(lmp->atom->mass[1], 4.0);
    ASSERT_DOUBLE_EQ(lmp->atom->mass[2], 2.4);
    ASSERT_EQ(lmp->atom->mass_setflag[1], 1);
    ASSERT_EQ(lmp->atom->mass_setflag[2], 1);
    ASSERT_EQ(lmp->atom->map_style, 1);
    ASSERT_EQ(lmp->atom->map_user, 1);
    ASSERT_EQ(lmp->atom->map_tag_max, 3);
    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("reset_ids");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(lmp->atom->map_tag_max, 2);
    ASSERT_EQ(lmp->atom->tag_consecutive(), 1);
}

TEST_F(AtomStyleTest, charge)
{
    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("atom_style charge");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_THAT(std::string(lmp->atom->atom_style), Eq("charge"));
    ASSERT_NE(lmp->atom->avec, nullptr);
    ASSERT_EQ(lmp->atom->natoms, 0);
    ASSERT_EQ(lmp->atom->nlocal, 0);
    ASSERT_EQ(lmp->atom->nghost, 0);
    ASSERT_EQ(lmp->atom->nmax, 1);
    ASSERT_EQ(lmp->atom->tag_enable, 1);
    ASSERT_EQ(lmp->atom->molecular, 0);
    ASSERT_EQ(lmp->atom->nellipsoids, 0);
    ASSERT_EQ(lmp->atom->nlines, 0);
    ASSERT_EQ(lmp->atom->ntris, 0);
    ASSERT_EQ(lmp->atom->nbodies, 0);
    ASSERT_EQ(lmp->atom->nbonds, 0);
    ASSERT_EQ(lmp->atom->nangles, 0);
    ASSERT_EQ(lmp->atom->ndihedrals, 0);
    ASSERT_EQ(lmp->atom->nimpropers, 0);
    ASSERT_EQ(lmp->atom->ntypes, 0);
    ASSERT_EQ(lmp->atom->nbondtypes, 0);
    ASSERT_EQ(lmp->atom->nangletypes, 0);
    ASSERT_EQ(lmp->atom->ndihedraltypes, 0);
    ASSERT_EQ(lmp->atom->nimpropertypes, 0);
    ASSERT_EQ(lmp->atom->bond_per_atom, 0);
    ASSERT_EQ(lmp->atom->angle_per_atom, 0);
    ASSERT_EQ(lmp->atom->dihedral_per_atom, 0);
    ASSERT_EQ(lmp->atom->improper_per_atom, 0);
    ASSERT_EQ(lmp->atom->extra_bond_per_atom, 0);
    ASSERT_EQ(lmp->atom->extra_angle_per_atom, 0);
    ASSERT_EQ(lmp->atom->extra_dihedral_per_atom, 0);
    ASSERT_EQ(lmp->atom->extra_improper_per_atom, 0);

    ASSERT_EQ(lmp->atom->sphere_flag, 0);
    ASSERT_EQ(lmp->atom->ellipsoid_flag, 0);
    ASSERT_EQ(lmp->atom->line_flag, 0);
    ASSERT_EQ(lmp->atom->tri_flag, 0);
    ASSERT_EQ(lmp->atom->body_flag, 0);
    ASSERT_EQ(lmp->atom->peri_flag, 0);
    ASSERT_EQ(lmp->atom->electron_flag, 0);
    ASSERT_EQ(lmp->atom->wavepacket_flag, 0);
    ASSERT_EQ(lmp->atom->sph_flag, 0);
    ASSERT_EQ(lmp->atom->molecule_flag, 0);
    ASSERT_EQ(lmp->atom->molindex_flag, 0);
    ASSERT_EQ(lmp->atom->molatom_flag, 0);
    ASSERT_EQ(lmp->atom->q_flag, 1);
    ASSERT_EQ(lmp->atom->mu_flag, 0);
    ASSERT_EQ(lmp->atom->rmass_flag, 0);
    ASSERT_EQ(lmp->atom->radius_flag, 0);
    ASSERT_EQ(lmp->atom->omega_flag, 0);
    ASSERT_EQ(lmp->atom->torque_flag, 0);
    ASSERT_EQ(lmp->atom->angmom_flag, 0);
    ASSERT_EQ(lmp->atom->vfrac_flag, 0);
    ASSERT_EQ(lmp->atom->spin_flag, 0);
    ASSERT_EQ(lmp->atom->eradius_flag, 0);
    ASSERT_EQ(lmp->atom->ervel_flag, 0);
    ASSERT_EQ(lmp->atom->erforce_flag, 0);
    ASSERT_EQ(lmp->atom->cs_flag, 0);
    ASSERT_EQ(lmp->atom->csforce_flag, 0);
    ASSERT_EQ(lmp->atom->vforce_flag, 0);
    ASSERT_EQ(lmp->atom->ervelforce_flag, 0);
    ASSERT_EQ(lmp->atom->etag_flag, 0);
    ASSERT_EQ(lmp->atom->rho_flag, 0);
    ASSERT_EQ(lmp->atom->esph_flag, 0);
    ASSERT_EQ(lmp->atom->cv_flag, 0);
    ASSERT_EQ(lmp->atom->vest_flag, 0);
    ASSERT_EQ(lmp->atom->dpd_flag, 0);
    ASSERT_EQ(lmp->atom->edpd_flag, 0);
    ASSERT_EQ(lmp->atom->tdpd_flag, 0);
    ASSERT_EQ(lmp->atom->mesont_flag, 0);
    ASSERT_EQ(lmp->atom->sp_flag, 0);
    ASSERT_EQ(lmp->atom->x0_flag, 0);
    ASSERT_EQ(lmp->atom->smd_flag, 0);
    ASSERT_EQ(lmp->atom->damage_flag, 0);
    ASSERT_EQ(lmp->atom->contact_radius_flag, 0);
    ASSERT_EQ(lmp->atom->smd_data_9_flag, 0);
    ASSERT_EQ(lmp->atom->smd_stress_flag, 0);
    ASSERT_EQ(lmp->atom->eff_plastic_strain_flag, 0);
    ASSERT_EQ(lmp->atom->eff_plastic_strain_rate_flag, 0);
    ASSERT_EQ(lmp->atom->pdscale, 1.0);

    ASSERT_NE(lmp->atom->tag, nullptr);
    ASSERT_NE(lmp->atom->type, nullptr);
    ASSERT_NE(lmp->atom->mask, nullptr);
    ASSERT_NE(lmp->atom->image, nullptr);
    ASSERT_NE(lmp->atom->x, nullptr);
    ASSERT_NE(lmp->atom->v, nullptr);
    ASSERT_NE(lmp->atom->f, nullptr);
    ASSERT_NE(lmp->atom->q, nullptr);
    ASSERT_EQ(lmp->atom->mu, nullptr);
    ASSERT_EQ(lmp->atom->omega, nullptr);
    ASSERT_EQ(lmp->atom->angmom, nullptr);
    ASSERT_EQ(lmp->atom->torque, nullptr);
    ASSERT_EQ(lmp->atom->radius, nullptr);
    ASSERT_EQ(lmp->atom->rmass, nullptr);
    ASSERT_EQ(lmp->atom->ellipsoid, nullptr);
    ASSERT_EQ(lmp->atom->line, nullptr);
    ASSERT_EQ(lmp->atom->tri, nullptr);
    ASSERT_EQ(lmp->atom->body, nullptr);
    ASSERT_EQ(lmp->atom->molecule, nullptr);
    ASSERT_EQ(lmp->atom->molindex, nullptr);
    ASSERT_EQ(lmp->atom->molatom, nullptr);
    ASSERT_EQ(lmp->atom->num_bond, nullptr);
    ASSERT_EQ(lmp->atom->bond_type, nullptr);
    ASSERT_EQ(lmp->atom->bond_atom, nullptr);
    ASSERT_EQ(lmp->atom->num_angle, nullptr);
    ASSERT_EQ(lmp->atom->angle_type, nullptr);
    ASSERT_EQ(lmp->atom->angle_atom1, nullptr);
    ASSERT_EQ(lmp->atom->angle_atom2, nullptr);
    ASSERT_EQ(lmp->atom->angle_atom3, nullptr);
    ASSERT_EQ(lmp->atom->num_dihedral, nullptr);
    ASSERT_EQ(lmp->atom->dihedral_type, nullptr);
    ASSERT_EQ(lmp->atom->dihedral_atom1, nullptr);
    ASSERT_EQ(lmp->atom->dihedral_atom2, nullptr);
    ASSERT_EQ(lmp->atom->dihedral_atom3, nullptr);
    ASSERT_EQ(lmp->atom->dihedral_atom4, nullptr);
    ASSERT_EQ(lmp->atom->num_improper, nullptr);
    ASSERT_EQ(lmp->atom->improper_type, nullptr);
    ASSERT_EQ(lmp->atom->improper_atom1, nullptr);
    ASSERT_EQ(lmp->atom->improper_atom2, nullptr);
    ASSERT_EQ(lmp->atom->improper_atom3, nullptr);
    ASSERT_EQ(lmp->atom->improper_atom4, nullptr);
    ASSERT_EQ(lmp->atom->maxspecial, 1);
    ASSERT_EQ(lmp->atom->nspecial, nullptr);
    ASSERT_EQ(lmp->atom->special, nullptr);
    ASSERT_EQ(lmp->atom->vfrac, nullptr);
    ASSERT_EQ(lmp->atom->s0, nullptr);
    ASSERT_EQ(lmp->atom->x0, nullptr);
    ASSERT_EQ(lmp->atom->sp, nullptr);
    ASSERT_EQ(lmp->atom->fm, nullptr);
    ASSERT_EQ(lmp->atom->fm_long, nullptr);
    ASSERT_EQ(lmp->atom->spin, nullptr);
    ASSERT_EQ(lmp->atom->eradius, nullptr);
    ASSERT_EQ(lmp->atom->ervel, nullptr);
    ASSERT_EQ(lmp->atom->erforce, nullptr);
    ASSERT_EQ(lmp->atom->ervelforce, nullptr);
    ASSERT_EQ(lmp->atom->cs, nullptr);
    ASSERT_EQ(lmp->atom->csforce, nullptr);
    ASSERT_EQ(lmp->atom->vforce, nullptr);
    ASSERT_EQ(lmp->atom->etag, nullptr);
    ASSERT_EQ(lmp->atom->uCond, nullptr);
    ASSERT_EQ(lmp->atom->uMech, nullptr);
    ASSERT_EQ(lmp->atom->uChem, nullptr);
    ASSERT_EQ(lmp->atom->uCG, nullptr);
    ASSERT_EQ(lmp->atom->uCGnew, nullptr);
    ASSERT_EQ(lmp->atom->duChem, nullptr);
    ASSERT_EQ(lmp->atom->dpdTheta, nullptr);
    ASSERT_EQ(lmp->atom->cc, nullptr);
    ASSERT_EQ(lmp->atom->cc_flux, nullptr);
    ASSERT_EQ(lmp->atom->edpd_temp, nullptr);
    ASSERT_EQ(lmp->atom->edpd_flux, nullptr);
    ASSERT_EQ(lmp->atom->edpd_cv, nullptr);
    ASSERT_EQ(lmp->atom->length, nullptr);
    ASSERT_EQ(lmp->atom->buckling, nullptr);
    ASSERT_EQ(lmp->atom->bond_nt, nullptr);
    ASSERT_EQ(lmp->atom->contact_radius, nullptr);
    ASSERT_EQ(lmp->atom->smd_data_9, nullptr);
    ASSERT_EQ(lmp->atom->smd_stress, nullptr);
    ASSERT_EQ(lmp->atom->eff_plastic_strain, nullptr);
    ASSERT_EQ(lmp->atom->eff_plastic_strain_rate, nullptr);
    ASSERT_EQ(lmp->atom->damage, nullptr);
    ASSERT_EQ(lmp->atom->rho, nullptr);
    ASSERT_EQ(lmp->atom->drho, nullptr);
    ASSERT_EQ(lmp->atom->esph, nullptr);
    ASSERT_EQ(lmp->atom->desph, nullptr);
    ASSERT_EQ(lmp->atom->cv, nullptr);
    ASSERT_EQ(lmp->atom->vest, nullptr);
    ASSERT_EQ(lmp->atom->nmolecule, 0);
    ASSERT_EQ(lmp->atom->molecules, nullptr);
    ASSERT_EQ(lmp->atom->nivector, 0);
    ASSERT_EQ(lmp->atom->ndvector, 0);
    ASSERT_EQ(lmp->atom->iname, nullptr);
    ASSERT_EQ(lmp->atom->dname, nullptr);
    ASSERT_EQ(lmp->atom->mass, nullptr);
    ASSERT_EQ(lmp->atom->mass_setflag, nullptr);
    ASSERT_EQ(lmp->atom->nextra_grow, 0);
    ASSERT_EQ(lmp->atom->nextra_restart, 0);
    ASSERT_EQ(lmp->atom->nextra_border, 0);
    ASSERT_EQ(lmp->atom->nextra_grow_max, 0);
    ASSERT_EQ(lmp->atom->nextra_restart_max, 0);
    ASSERT_EQ(lmp->atom->nextra_border_max, 0);
    ASSERT_EQ(lmp->atom->nextra_store, 0);
    ASSERT_EQ(lmp->atom->extra_grow, nullptr);
    ASSERT_EQ(lmp->atom->extra_restart, nullptr);
    ASSERT_EQ(lmp->atom->extra_border, nullptr);
    ASSERT_EQ(lmp->atom->extra, nullptr);
    ASSERT_EQ(lmp->atom->sametag, nullptr);
    ASSERT_EQ(lmp->atom->map_style, 0);
    ASSERT_EQ(lmp->atom->map_user, 0);
    ASSERT_EQ(lmp->atom->map_tag_max, -1);

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("create_box 2 box");
    lmp->input->one("create_atoms 1 single -2.0  2.0  0.1");
    lmp->input->one("create_atoms 1 single -2.0 -2.0 -0.1");
    lmp->input->one("create_atoms 2 single  2.0  2.0 -0.1");
    lmp->input->one("create_atoms 2 single  2.0 -2.0  0.1");
    lmp->input->one("mass 1 4.0");
    lmp->input->one("mass 2 2.4");
    lmp->input->one("set atom 1 charge -0.5");
    lmp->input->one("set atom 2 charge  0.5");
    lmp->input->one("set atom 3 charge -1.0");
    lmp->input->one("set atom 4 charge  1.0");
    lmp->input->one("pair_coeff * *");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_THAT(std::string(lmp->atom->atom_style), Eq("charge"));
    ASSERT_NE(lmp->atom->avec, nullptr);
    ASSERT_EQ(lmp->atom->natoms, 4);
    ASSERT_EQ(lmp->atom->nlocal, 4);
    ASSERT_EQ(lmp->atom->nghost, 0);
    ASSERT_NE(lmp->atom->nmax, -1);
    ASSERT_EQ(lmp->atom->tag_enable, 1);
    ASSERT_EQ(lmp->atom->molecular, 0);
    ASSERT_EQ(lmp->atom->ntypes, 2);

    ASSERT_NE(lmp->atom->mass, nullptr);
    ASSERT_NE(lmp->atom->mass_setflag, nullptr);
    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("pair_coeff * *");
    lmp->input->one("write_data test_atom_styles.data nocoeff");
    lmp->input->one("clear");
    lmp->input->one("atom_style charge");
    lmp->input->one("pair_style zero 4.0");
    lmp->input->one("units real");
    lmp->input->one("atom_modify map array");
    lmp->input->one("read_data test_atom_styles.data");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_THAT(std::string(lmp->atom->atom_style), Eq("charge"));
    ASSERT_NE(lmp->atom->avec, nullptr);
    ASSERT_EQ(lmp->atom->natoms, 4);
    ASSERT_EQ(lmp->atom->nlocal, 4);
    ASSERT_EQ(lmp->atom->nghost, 0);
    ASSERT_NE(lmp->atom->nmax, -1);
    ASSERT_EQ(lmp->atom->tag_enable, 1);
    ASSERT_EQ(lmp->atom->molecular, 0);
    ASSERT_EQ(lmp->atom->ntypes, 2);
    ASSERT_EQ(lmp->atom->q_flag, 1);
    ASSERT_NE(lmp->atom->sametag, nullptr);
    ASSERT_EQ(lmp->atom->tag_consecutive(), 1);
    ASSERT_EQ(lmp->atom->map_style, 1);
    ASSERT_EQ(lmp->atom->map_user, 1);
    ASSERT_EQ(lmp->atom->map_tag_max, 4);

    double **x = lmp->atom->x;
    double **v = lmp->atom->v;
    double *q  = lmp->atom->q;
    ASSERT_DOUBLE_EQ(x[0][0], -2.0);
    ASSERT_DOUBLE_EQ(x[0][1], 2.0);
    ASSERT_DOUBLE_EQ(x[0][2], 0.1);
    ASSERT_DOUBLE_EQ(x[1][0], -2.0);
    ASSERT_DOUBLE_EQ(x[1][1], -2.0);
    ASSERT_DOUBLE_EQ(x[1][2], -0.1);
    ASSERT_DOUBLE_EQ(x[2][0], 2.0);
    ASSERT_DOUBLE_EQ(x[2][1], 2.0);
    ASSERT_DOUBLE_EQ(x[2][2], -0.1);
    ASSERT_DOUBLE_EQ(x[3][0], 2.0);
    ASSERT_DOUBLE_EQ(x[3][1], -2.0);
    ASSERT_DOUBLE_EQ(x[3][2], 0.1);
    ASSERT_DOUBLE_EQ(v[0][0], 0.0);
    ASSERT_DOUBLE_EQ(v[0][1], 0.0);
    ASSERT_DOUBLE_EQ(v[0][2], 0.0);
    ASSERT_DOUBLE_EQ(v[1][0], 0.0);
    ASSERT_DOUBLE_EQ(v[1][1], 0.0);
    ASSERT_DOUBLE_EQ(v[1][2], 0.0);
    ASSERT_DOUBLE_EQ(v[2][0], 0.0);
    ASSERT_DOUBLE_EQ(v[2][1], 0.0);
    ASSERT_DOUBLE_EQ(v[2][2], 0.0);
    ASSERT_DOUBLE_EQ(v[3][0], 0.0);
    ASSERT_DOUBLE_EQ(v[3][1], 0.0);
    ASSERT_DOUBLE_EQ(v[3][2], 0.0);
    ASSERT_DOUBLE_EQ(q[0], -0.5);
    ASSERT_DOUBLE_EQ(q[1], 0.5);
    ASSERT_DOUBLE_EQ(q[2], -1.0);
    ASSERT_DOUBLE_EQ(q[3],  1.0);

    ASSERT_DOUBLE_EQ(lmp->atom->mass[1], 4.0);
    ASSERT_DOUBLE_EQ(lmp->atom->mass[2], 2.4);
    ASSERT_EQ(lmp->atom->mass_setflag[1], 1);
    ASSERT_EQ(lmp->atom->mass_setflag[2], 1);

    if (!verbose) ::testing::internal::CaptureStdout();
    lmp->input->one("pair_coeff * *");
    lmp->input->one("group two id 2:4:2");
    lmp->input->one("delete_atoms group two compress no");
    lmp->input->one("write_restart test_atom_styles.restart");
    lmp->input->one("clear");
    lmp->input->one("read_restart test_atom_styles.restart");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_THAT(std::string(lmp->atom->atom_style), Eq("charge"));
    ASSERT_NE(lmp->atom->avec, nullptr);
    ASSERT_EQ(lmp->atom->natoms, 2);
    ASSERT_EQ(lmp->atom->nlocal, 2);
    ASSERT_EQ(lmp->atom->nghost, 0);
    ASSERT_NE(lmp->atom->nmax, -1);
    ASSERT_EQ(lmp->atom->tag_enable, 1);
    ASSERT_EQ(lmp->atom->molecular, 0);
    ASSERT_EQ(lmp->atom->ntypes, 2);
    ASSERT_EQ(lmp->atom->tag_consecutive(), 0);
    ASSERT_EQ(lmp->atom->map_tag_max, 3);

    x = lmp->atom->x;
    v = lmp->atom->v;
    q = lmp->atom->q;
    ASSERT_DOUBLE_EQ(x[0][0], -2.0);
    ASSERT_DOUBLE_EQ(x[0][1], 2.0);
    ASSERT_DOUBLE_EQ(x[0][2], 0.1);
    ASSERT_DOUBLE_EQ(x[1][0], 2.0);
    ASSERT_DOUBLE_EQ(x[1][1], 2.0);
    ASSERT_DOUBLE_EQ(x[1][2], -0.1);
    ASSERT_DOUBLE_EQ(v[0][0], 0.0);
    ASSERT_DOUBLE_EQ(v[0][1], 0.0);
    ASSERT_DOUBLE_EQ(v[0][2], 0.0);
    ASSERT_DOUBLE_EQ(v[1][0], 0.0);
    ASSERT_DOUBLE_EQ(v[1][1], 0.0);
    ASSERT_DOUBLE_EQ(v[1][2], 0.0);
    ASSERT_DOUBLE_EQ(q[0], -0.5);
    ASSERT_DOUBLE_EQ(q[1], -1.0);

    ASSERT_DOUBLE_EQ(lmp->atom->mass[1], 4.0);
    ASSERT_DOUBLE_EQ(lmp->atom->mass[2], 2.4);
    ASSERT_EQ(lmp->atom->mass_setflag[1], 1);
    ASSERT_EQ(lmp->atom->mass_setflag[2], 1);
}

} // namespace LAMMPS_NS

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleMock(&argc, argv);
    if ((argc > 1) && (strcmp(argv[1], "-v") == 0)) verbose = true;

    int rv = RUN_ALL_TESTS();
    MPI_Finalize();
    return rv;
}
