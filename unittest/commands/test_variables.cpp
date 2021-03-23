/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://lammps.sandia.gov/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "lammps.h"

#include "atom.h"
#include "domain.h"
#include "group.h"
#include "info.h"
#include "input.h"
#include "math_const.h"
#include "region.h"
#include "variable.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <cstring>
#include <vector>

// whether to print verbose output (i.e. not capturing LAMMPS screen output).
bool verbose = false;

#if defined(OMPI_MAJOR_VERSION)
const bool have_openmpi = true;
#else
const bool have_openmpi = false;
#endif

using LAMMPS_NS::MathConst::MY_PI;
using LAMMPS_NS::utils::split_words;

namespace LAMMPS_NS {
using ::testing::ExitedWithCode;
using ::testing::MatchesRegex;
using ::testing::StrEq;

#define TEST_FAILURE(errmsg, ...)                                 \
    if (Info::has_exceptions()) {                                 \
        ::testing::internal::CaptureStdout();                     \
        ASSERT_ANY_THROW({__VA_ARGS__});                          \
        auto mesg = ::testing::internal::GetCapturedStdout();     \
        if (verbose) std::cout << mesg;                           \
        ASSERT_THAT(mesg, MatchesRegex(errmsg));                  \
    } else {                                                      \
        if (!have_openmpi) {                                      \
            ::testing::internal::CaptureStdout();                 \
            ASSERT_DEATH({__VA_ARGS__}, "");                      \
            auto mesg = ::testing::internal::GetCapturedStdout(); \
            if (verbose) std::cout << mesg;                       \
            ASSERT_THAT(mesg, MatchesRegex(errmsg));              \
        }                                                         \
    }

class VariableTest : public ::testing::Test {
protected:
    LAMMPS *lmp;
    Group *group;
    Domain *domain;
    Variable *variable;

    void SetUp() override
    {
        const char *args[] = {"VariableTest", "-log", "none", "-echo", "screen", "-nocite"};
        char **argv        = (char **)args;
        int argc           = sizeof(args) / sizeof(char *);
        if (!verbose) ::testing::internal::CaptureStdout();
        lmp = new LAMMPS(argc, argv, MPI_COMM_WORLD);
        if (!verbose) ::testing::internal::GetCapturedStdout();
        group    = lmp->group;
        domain   = lmp->domain;
        variable = lmp->input->variable;
    }

    void TearDown() override
    {
        if (!verbose) ::testing::internal::CaptureStdout();
        delete lmp;
        if (!verbose) ::testing::internal::GetCapturedStdout();
        std::cout.flush();
        unlink("test_variable.file");
        unlink("test_variable.atomfile");
    }

    void command(const std::string &cmd) { lmp->input->one(cmd); }

    void atomic_system()
    {
        if (!verbose) ::testing::internal::CaptureStdout();
        command("units real");
        command("lattice sc 1.0 origin 0.125 0.125 0.125");
        command("region box block -2 2 -2 2 -2 2");
        command("create_box 8 box");
        command("create_atoms 1 box");
        command("mass * 1.0");
        command("region left block -2.0 -1.0 INF INF INF INF");
        command("region right block 0.5  2.0 INF INF INF INF");
        command("region top block INF INF -2.0 -1.0 INF INF");
        command("set region left type 2");
        command("set region right type 3");
        if (!verbose) ::testing::internal::GetCapturedStdout();
    }

    void molecular_system()
    {
        if (!verbose) ::testing::internal::CaptureStdout();
        command("fix props all property/atom mol rmass q");
        if (!verbose) ::testing::internal::GetCapturedStdout();
        atomic_system();
        if (!verbose) ::testing::internal::CaptureStdout();
        command("variable molid atom floor(id/4)+1");
        command("variable charge atom 2.0*sin(PI/32*id)");
        command("set atom * mol v_molid");
        command("set atom * charge v_charge");
        command("set type 1 mass 0.5");
        command("set type 2*4 mass 2.0");
        if (!verbose) ::testing::internal::GetCapturedStdout();
    }

    void file_vars()
    {
        FILE *fp = fopen("test_variable.file", "w");
        fputs("# test file for file style variable\n\n\none\n  two  \n\n"
              "three  # with comment\nfour   ! with non-comment\n"
              "# comments only\n	five\n#END\n",
              fp);
        fclose(fp);
        fp = fopen("test_variable.atomfile", "w");

        fputs("# test file for atomfile style variable\n\n"
              "4  # four lines\n4 0.5   #with comment\n"
              "2 -0.5         \n3 1.5\n1 -1.5\n\n"
              "2\n10 1.0 # test\n13 1.0\n\n######\n"
              "4\n1 4.0 # test\n2 3.0\n3 2.0\n4 1.0\n#END\n",
              fp);
        fclose(fp);
    }
};

TEST_F(VariableTest, CreateDelete)
{
    file_vars();
    ASSERT_EQ(variable->nvar, 0);
    if (!verbose) ::testing::internal::CaptureStdout();
    command("variable one    index     1 2 3 4");
    command("variable two    equal     1");
    command("variable two    equal     2");
    command("variable three  string    four");
    command("variable three  string    three");
    command("variable four1  loop      4");
    command("variable four2  loop      2 4");
    command("variable five1  loop      100 pad");
    command("variable five2  loop      100 200 pad");
    command("variable six    world     one");
    command("variable seven  format    two \"%5.2f\"");
    command("variable eight  getenv    PWD");
    command("variable eight  getenv    HOME");
    command("variable nine   file      test_variable.file");
    command("variable dummy  index     0");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(variable->nvar, 12);
    if (!verbose) ::testing::internal::CaptureStdout();
    command("variable dummy  delete");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(variable->nvar, 11);

    TEST_FAILURE(".*ERROR: Illegal variable command.*", command("variable"););
    TEST_FAILURE(".*ERROR: Illegal variable command.*", command("variable dummy index"););
    TEST_FAILURE(".*ERROR: Illegal variable command.*", command("variable dummy delete xxx"););
    TEST_FAILURE(".*ERROR: Cannot redefine variable as a different style.*",
                 command("variable two string xxx"););
    TEST_FAILURE(".*ERROR: Cannot redefine variable as a different style.*",
                 command("variable two getenv xxx"););
    TEST_FAILURE(".*ERROR: Cannot redefine variable as a different style.*",
                 command("variable one equal 2"););
    TEST_FAILURE(".*ERROR: Cannot use atomfile-style variable unless an atom map exists.*",
                 command("variable ten    atomfile  test_variable.atomfile"););
    TEST_FAILURE(".*ERROR on proc 0: Cannot open file variable file test_variable.xxx.*",
                 command("variable nine1  file      test_variable.xxx"););
}

TEST_F(VariableTest, AtomicSystem)
{
    command("atom_modify map array");
    atomic_system();
    file_vars();

    if (!verbose) ::testing::internal::CaptureStdout();
    command("variable  one  index     1 2 3 4");
    command("variable  id   atom      type");
    command("variable  id   atom      id");
    command("variable  ten  atomfile  test_variable.atomfile");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(variable->nvar, 3);

    TEST_FAILURE(".*ERROR: Cannot redefine variable as a different style.*",
                 command("variable one atom x"););
    TEST_FAILURE(".*ERROR on proc 0: Cannot open file variable file test_variable.xxx.*",
                 command("variable ten1   atomfile  test_variable.xxx"););
}

TEST_F(VariableTest, Expressions)
{
    atomic_system();
    ASSERT_EQ(variable->nvar, 0);
    if (!verbose) ::testing::internal::CaptureStdout();
    command("variable one    index     1");
    command("variable two    equal     2");
    command("variable three  equal     v_one+v_two");
    command("variable four   equal     PI");
    command("variable five   equal     version");
    command("variable six    equal     XXX");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(variable->nvar, 6);

    int ivar = variable->find("one");
    ASSERT_FALSE(variable->equalstyle(ivar));
    ivar = variable->find("two");
    ASSERT_TRUE(variable->equalstyle(ivar));
    ASSERT_DOUBLE_EQ(variable->compute_equal(ivar), 2.0);
    ivar = variable->find("three");
    ASSERT_DOUBLE_EQ(variable->compute_equal(ivar), 3.0);
    ivar = variable->find("four");
    ASSERT_DOUBLE_EQ(variable->compute_equal(ivar), MY_PI);
    ivar = variable->find("five");
    ASSERT_GE(variable->compute_equal(ivar), 20210310);

    TEST_FAILURE(".*ERROR: Variable six: Invalid thermo keyword 'XXX' in variable formula.*",
                 command("print \"${six}\""););
}

TEST_F(VariableTest, Functions)
{
    atomic_system();
    file_vars();

    ASSERT_EQ(variable->nvar, 0);
    if (!verbose) ::testing::internal::CaptureStdout();
    command("variable one    index     1");
    command("variable two    equal     random(1,2,643532)");
    command("variable three  equal     atan2(v_one,1)");
    command("variable four   equal     atan2()");
    if (!verbose) ::testing::internal::GetCapturedStdout();
    ASSERT_EQ(variable->nvar, 4);

    int ivar = variable->find("two");
    ASSERT_GT(variable->compute_equal(ivar), 0.99);
    ASSERT_LT(variable->compute_equal(ivar), 2.01);
    ivar = variable->find("three");
    ASSERT_DOUBLE_EQ(variable->compute_equal(ivar), 0.25 * MY_PI);
    TEST_FAILURE(".*ERROR: Variable four: Invalid syntax in variable formula.*",
                 command("print \"${four}\""););
}
} // namespace LAMMPS_NS

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleMock(&argc, argv);

    if (have_openmpi && !LAMMPS_NS::Info::has_exceptions())
        std::cout << "Warning: using OpenMPI without exceptions. "
                     "Death tests will be skipped\n";

    // handle arguments passed via environment variable
    if (const char *var = getenv("TEST_ARGS")) {
        std::vector<std::string> env = split_words(var);
        for (auto arg : env) {
            if (arg == "-v") {
                verbose = true;
            }
        }
    }

    if ((argc > 1) && (strcmp(argv[1], "-v") == 0)) verbose = true;

    int rv = RUN_ALL_TESTS();
    MPI_Finalize();
    return rv;
}
