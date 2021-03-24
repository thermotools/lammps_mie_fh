
from __future__ import print_function

def square(val):
    return val*val

def bool_to_val(txt):
    if txt.upper() in ["TRUE", "YES"]:
        return 1.0
    return 0.0

def printnum():
    print("2.25")

def printtxt():
    print("sometext")

def getidxvar(lmpptr):
    from lammps import lammps
    lmp = lammps(ptr=lmpptr)
    val = lmp.extract_variable("idx")
    print(val)
