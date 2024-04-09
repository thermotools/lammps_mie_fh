# Feynman-Hibbs Corrected Mie Pair Potential for LAMMPS

This repository contains an extension for the Large-scale Atomic/Molecular Massively Parallel Simulator (LAMMPS) that implements a quantum correction pair potential based on Feynman-Hibbs corrected Mie potential.

Table of Contents
-----------------

* [Introduction](#introduction)
* [Installation](#installation)
* [Usage](#usage)
* [Citation](#citation)

Introduction
------------
### Mie-FH Potential

The Mie-FH potential is defined as:

[![\\ \begin{aligned} \\ \frac{u_{i j}\left(r_{i j}\right)}{\mathcal{C}\left(\gamma_{r, i j}, \gamma_{a, i j}\right) \epsilon_{i j}}= & \frac{\sigma_{i j}^{\gamma_r, i j}}{r_{i j}^{\gamma_r, i j}}-\frac{\sigma_{i j}^{\gamma_a, i j}}{r_{i j}^{\gamma_{a, i j}}} \\ \\ & +D\left(Q_1\left(\gamma_{r, i j}\right) \frac{\sigma_{i j}^{\gamma_r, i j+2}}{r_{i j}^{\gamma_r, i j}+2}-Q_1\left(\gamma_{a, i j}\right) \frac{\sigma_{i j}^{\gamma_{a, i j+2}}}{r_{i j}^{\gamma_{a, i j}+2}}\right) \\ \\ & +D^2\left(Q_2\left(\gamma_{r, i j}\right) \frac{\sigma_{i j}^{\gamma_r, i j+4}}{r_{i j}^{\gamma_r, i j}+4}-Q_2\left(\gamma_{a, i j}\right) \frac{\sigma_{i j}^{\gamma_a, i j+4}}{r_{i j}^{\gamma_{a, i j}+4}}\right) \\ \end{aligned}](https://latex.codecogs.com/svg.latex?%5C%5C%20%5Cbegin%7Baligned%7D%20%5C%5C%20%5Cfrac%7Bu_%7Bi%20j%7D%5Cleft(r_%7Bi%20j%7D%5Cright)%7D%7B%5Cmathcal%7BC%7D%5Cleft(%5Cgamma_%7Br%2C%20i%20j%7D%2C%20%5Cgamma_%7Ba%2C%20i%20j%7D%5Cright)%20%5Cepsilon_%7Bi%20j%7D%7D%3D%20%26%20%5Cfrac%7B%5Csigma_%7Bi%20j%7D%5E%7B%5Cgamma_r%2C%20i%20j%7D%7D%7Br_%7Bi%20j%7D%5E%7B%5Cgamma_r%2C%20i%20j%7D%7D-%5Cfrac%7B%5Csigma_%7Bi%20j%7D%5E%7B%5Cgamma_a%2C%20i%20j%7D%7D%7Br_%7Bi%20j%7D%5E%7B%5Cgamma_%7Ba%2C%20i%20j%7D%7D%7D%20%5C%5C%20%5C%5C%20%26%20%2BD%5Cleft(Q_1%5Cleft(%5Cgamma_%7Br%2C%20i%20j%7D%5Cright)%20%5Cfrac%7B%5Csigma_%7Bi%20j%7D%5E%7B%5Cgamma_r%2C%20i%20j%2B2%7D%7D%7Br_%7Bi%20j%7D%5E%7B%5Cgamma_r%2C%20i%20j%7D%2B2%7D-Q_1%5Cleft(%5Cgamma_%7Ba%2C%20i%20j%7D%5Cright)%20%5Cfrac%7B%5Csigma_%7Bi%20j%7D%5E%7B%5Cgamma_%7Ba%2C%20i%20j%2B2%7D%7D%7D%7Br_%7Bi%20j%7D%5E%7B%5Cgamma_%7Ba%2C%20i%20j%7D%2B2%7D%7D%5Cright)%20%5C%5C%20%5C%5C%20%26%20%2BD%5E2%5Cleft(Q_2%5Cleft(%5Cgamma_%7Br%2C%20i%20j%7D%5Cright)%20%5Cfrac%7B%5Csigma_%7Bi%20j%7D%5E%7B%5Cgamma_r%2C%20i%20j%2B4%7D%7D%7Br_%7Bi%20j%7D%5E%7B%5Cgamma_r%2C%20i%20j%7D%2B4%7D-Q_2%5Cleft(%5Cgamma_%7Ba%2C%20i%20j%7D%5Cright)%20%5Cfrac%7B%5Csigma_%7Bi%20j%7D%5E%7B%5Cgamma_a%2C%20i%20j%2B4%7D%7D%7Br_%7Bi%20j%7D%5E%7B%5Cgamma_%7Ba%2C%20i%20j%7D%2B4%7D%7D%5Cright)%20%5C%5C%20%5Cend%7Baligned%7D)](#_)

Installation
------------
The two new pair style files (mie_fh1/cut and mie_fh2/cut) are located at EXTRA-PAIR folder.
This extension should work with different versions of LAMMPS. However, we include the LAMMPS stable 2022 version here.

To install the Mie-FH, follow these steps:
```
git clone https://github.com/thermotools/lammps_mie_fh.git
cd lammps_mie_fh
cd src
make yes-extra-pair
make yes-molecule
make -j4 mpi
mpirun -np 4 lmp_mpi -in ../Mie-FH1-npt.lmp
```

Usage
-----

To use this potential, see example Mie-FH1-npt.lmp file. Note that the the Mie-FH pair potential is temperature dependant. 

Citation
--------

Please cite the following papers when using this extension in your research:

(**Mie**) Mie, G. (1903). Zur kinetischen Theorie der einatomigen Körper. Annalen der Physik, 316(8), 657-697.

(**Avendano**) Avendano, C., Lafitte, T., Galindo, A., Adjiman, C. S., Jackson, G., & Müller, E. A. (2011). SAFT-γ force field for the simulation of molecular fluids. 1. A single-site coarse grained model of carbon dioxide. The Journal of Physical Chemistry B, 115(38), 11154-11169.

(**Aasen**) Aasen, A., Hammer, M., Ervik, Å., Müller, E. A., & Wilhelmsen, Ø. (2019). Equation of state and force fields for Feynman–Hibbs-corrected Mie fluids. I. Application to pure helium, neon, hydrogen, and deuterium. The Journal of Chemical Physics, 151(6).

(**Trinh**) Trinh, Thuat T and Hammer, Morten and Sharma, Vishist and Wilhelmsen, Øivind (2024). Mie-FH: A quantum corrected pair potential in the LAMMPS simulation package for hydrogen mixtures. SoftwareX, 26, 101716.
```{bibtex}
@article{trinh2024mie,
  title={Mie--FH: A quantum corrected pair potential in the LAMMPS simulation package for hydrogen mixtures},
  author={Trinh, Thuat T and Hammer, Morten and Sharma, Vishist and Wilhelmsen, {\O}ivind},
  journal={SoftwareX},
  volume={26},
  pages={101716},
  year={2024},
  publisher={Elsevier}
}
```


If you have any questions or issues with this extension, please feel free to open an issue or contact the author directly.

**Happy simulations!**

