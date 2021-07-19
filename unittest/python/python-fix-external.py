import sys,os,unittest
from ctypes import *
from lammps import lammps, LMP_STYLE_GLOBAL, LMP_TYPE_VECTOR

# add timestep dependent force
def callback_one(lmp, ntimestep, nlocal, tag, x, f):
    lmp.fix_external_set_virial_global("ext",[1.0, 1.0, 1.0, 0.0, 0.0, 0.0])
    for i in range(nlocal):
        f[i][0] = float(ntimestep)
        f[i][1] = float(ntimestep)
        f[i][2] = float(ntimestep)
    if ntimestep < 10:
        lmp.fix_external_set_energy_global("ext", 0.5)
        lmp.fix_external_set_vector("ext", 1, ntimestep)
        lmp.fix_external_set_vector("ext", 3, 1.0)
        lmp.fix_external_set_vector("ext", 4, -0.25)
    else:
        lmp.fix_external_set_energy_global("ext", 1.0)
        lmp.fix_external_set_vector("ext", 2, ntimestep)
        lmp.fix_external_set_vector("ext", 5, -1.0)
        lmp.fix_external_set_vector("ext", 6, 0.25)

class PythonExternal(unittest.TestCase):
    def testExternalCallback(self):
        """Test fix external from Python with pf/callback"""

        machine=None
        if 'LAMMPS_MACHINE_NAME' in os.environ:
            machine=os.environ['LAMMPS_MACHINE_NAME']
        lmp=lammps(name=machine, cmdargs=['-nocite', '-log','none', '-echo', 'screen'])

        # a few commands to set up simple system
        basic_system="""lattice sc 1.0
                        region box block -1 1 -1 1 -1 1
                        create_box 1 box
                        create_atoms 1 box
                        mass 1 1.0
                        pair_style zero 0.1
                        pair_coeff 1 1
                        velocity all set 0.1 0.0 -0.1
                        thermo_style custom step temp pe ke etotal press
                        thermo 5
                        fix 1 all nve
                        fix ext all external pf/callback 5 1
                        fix_modify ext energy yes virial yes
"""
        lmp.commands_string(basic_system)
        lmp.fix_external_set_vector_length("ext",6);
        lmp.set_fix_external_callback("ext",callback_one,lmp)
        lmp.command("run 10 post no")
        self.assertAlmostEqual(lmp.get_thermo("temp"),1.0/30.0,14)
        self.assertAlmostEqual(lmp.get_thermo("pe"),1.0/8.0,14)
        self.assertAlmostEqual(lmp.get_thermo("press"),0.15416666666666667,14)
        val = 0.0
        for i in range(0,6):
            val += lmp.extract_fix("ext",LMP_STYLE_GLOBAL,LMP_TYPE_VECTOR,nrow=i)
        self.assertAlmostEqual(val,15.0,14)

    def testExternalArray(self):
        """Test fix external from Python with pf/array"""

        machine=None
        if 'LAMMPS_MACHINE_NAME' in os.environ:
            machine=os.environ['LAMMPS_MACHINE_NAME']
        lmp=lammps(name=machine, cmdargs=['-nocite', '-log','none', '-echo', 'screen'])

        # a few commands to set up simple system
        basic_system="""lattice sc 1.0
                        region box block -1 1 -1 1 -1 1
                        create_box 1 box
                        create_atoms 1 box
                        mass 1 1.0
                        pair_style zero 0.1
                        pair_coeff 1 1
                        velocity all set 0.1 0.0 -0.1
                        fix 1 all nve
                        fix ext all external pf/array 1
                        thermo_style custom step temp pe ke
                        thermo 5
"""
        lmp.commands_string(basic_system)
        force = lmp.fix_external_get_force("ext");
        nlocal = lmp.extract_setting("nlocal");
        for i in range(nlocal):
            force[i][0] = 0.0
            force[i][1] = 0.0
            force[i][2] = 0.0
        lmp.fix_external_set_energy_global("ext", 0.5)
            
        lmp.command("run 5 post no")
        self.assertAlmostEqual(lmp.get_thermo("temp"),4.0/525.0,14)
        self.assertAlmostEqual(lmp.get_thermo("pe"),1.0/16.0,14)

        force = lmp.fix_external_get_force("ext");
        nlocal = lmp.extract_setting("nlocal");
        for i in range(nlocal):
            force[i][0] = 6.0
            force[i][1] = 6.0
            force[i][2] = 6.0
        lmp.fix_external_set_energy_global("ext", 1.0)
        lmp.command("run 5 post no")
        self.assertAlmostEqual(lmp.get_thermo("temp"),1.0/30.0,14)
        self.assertAlmostEqual(lmp.get_thermo("pe"),1.0/8.0,14)

##############################
if __name__ == "__main__":
    unittest.main()

