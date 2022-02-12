#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Python version of the code Zhou04_create_v2.f
original author: X. W. Zhou, xzhou@sandia.gov
based on updates by: Lucas Hale lucas.hale@nist.gov
written by: Germain Clavier g.m.g.c.clavier@tue.nl
This script requires the numpy library.
"""

import sys
import argparse as ap
import numpy as np
from eamDatabase import Database
from datetime import date

def prof(at, r):
    atom = Database[at]
    f = np.zeros(r.shape)
    f = atom.fe * np.exp(-atom.beta1 * (r[r >= 0.5] / atom.re - 1.0))
    f = f / (1.0 + (r / atom.re - atom.ramda1) ** 20)
    return f


def pair(at1, at2, r):
    atom = Database[at1]
    psi1 = atom.A * np.exp(-atom.alpha * (r / atom.re - 1.0))
    psi1 /= 1.0 + (r / atom.re - atom.cai) ** 20
    psi2 = atom.B * np.exp(-atom.beta * (r / atom.re - 1.0))
    psi2 /= 1.0 + (r / atom.re - atom.ramda) ** 20
    if at1 == at2:
        psi = psi1 - psi2
    else:
        psia = psi1 - psi2
        atom2 = Database[at2]
        psi1 = atom2.A * np.exp(-atom2.alpha * (r / atom2.re - 1.0))
        psi1 /= 1.0 + (r / atom2.re - atom2.cai) ** 20
        psi2 = atom2.B * np.exp(-atom2.beta * (r / atom2.re - 1.0))
        psi2 /= 1.0 + (r / atom2.re - atom2.ramda) ** 20
        psib = psi1 - psi2
        prof1 = prof(at1, r)
        prof2 = prof(at2, r)
        psi = 0.5 * (prof2 / prof1 * psia + prof1 / prof2 * psib)
    return psi


def embed(at, rho):
    atom = Database[at]
    Fm33 = np.zeros(rho.shape)
    Fm33[rho < atom.rhoe] = atom.Fm3
    Fm33[rho >= atom.rhoe] = atom.Fm4
    emb = np.zeros(rho.shape)
    for i, r in enumerate(rho):
        if r == 0:
            emb[i] = 0
        elif r < atom.rhoin:
            dr = r / atom.rhoin - 1
            emb[i] = atom.Fi0 + atom.Fi1 * dr + atom.Fi2 * dr**2 + atom.Fi3 * dr**3
        elif r < atom.rhoout:
            dr = r / atom.rhoe - 1
            emb[i] = atom.Fm0 + atom.Fm1 * dr + atom.Fm2 * dr**2 + Fm33[i] * dr**3
        else:
            dr = r / atom.rhos
            emb[i] = atom.Fn * (1.0 - atom.fnn * np.log(dr)) * dr**atom.fnn
    return emb


def write_file(attypes, filename, Fr, rhor, z2r, nrho, drho, nr, dr, rc):
    struc = "fcc"
    with open(filename, "w") as f:
        f.write("DATE: " + date.today().strftime("%Y-%m-%d"))
        f.write(" CONTRIBUTORS: Xiaowang Zhou xzhou@sandia.gov and")
        f.write(" Lucas Hale lucas.hale@nist.gov")
        f.write(" Germain Clavier g.m.g.c.clavier@tue.nl/germain.clavier@gmail.com\n")
        f.write(" CITATION: X. W. Zhou, R. A. Johnson, H. N. G. Wadley, Phys. Rev. B, 69, 144113(2004)\n")
        f.write("Generated from create_eam.py, Python version of Zhou04_create_v2.f\n")
        f.write("{:<5d} {:<24}\n".format(len(attypes), " ".join(attypes)))
        f.write(
            "{:<5d} {:<24.16e} {:<5d} {:<24.16e} {:<24.16e}\n".format(
                nrho, drho, nr, dr, rc
            )
        )
        for at in attypes:
            atom = Database[at]
            f.write(
                "{:>5d} {:>15.5f} {:>15.5f} {:>8}\n".format(
                    atom.ielement, atom.amass, atom.blat, struc
                )
            )
            for i, fr in enumerate(Fr[at]):
                f.write(" {:>24.16E}".format(fr))
                if not (i + 1) % 5:
                    f.write("\n")
            for i, rho in enumerate(rhor[at]):
                f.write(" {:>24.16E}".format(rho))
                if not (i + 1) % 5:
                    f.write("\n")
        for n1 in range(len(attypes)):
            for n2 in range(n1 + 1):
                for i, z in enumerate(z2r[n1, n2]):
                    f.write(" {:>24.16E}".format(z))
                    if not (i + 1) % 5:
                        f.write("\n")

    return


def main():

    parser = ap.ArgumentParser(description="Script to make EAM alloy file inputs.")
    parser.add_argument("-n", "--names", dest="names", nargs="+", help="Atom names.")
    parser.add_argument("-nr", dest="nr", type=int, default=2000, help="Number of point in r space [default 2000].")
    parser.add_argument("-nrho", dest="nrho", type=int, default=2000, help="Number of point in rho space [default 2000].")
    args = parser.parse_args()
    atnames = args.names
    nr = args.nr
    nrho = args.nrho

    for n in atnames:
        try:
            Database[n]
        except KeyError:
            output = "Atom {} not found.\n".format(n)
            valid = "Valid inputs are: {}".format(" ".join(Database.keys()))
            sys.exit("".join([output, valid]))

    ntypes = len(atnames)
    outfilename = "".join([*atnames, ".set"])
    rhor = {}
    Fr = {}

    alatmax = max([Database[at].blat for at in atnames])
    rhoemax = max([Database[at].rhoe for at in atnames])
    rc = np.sqrt(10.0) / 2.0 * alatmax
    rst = 0.5
    r = np.linspace(0.0, rc, num=nr, dtype=np.double)
    dr = r[1] - r[0]
    r[r < rst] = rst
    z2r = np.zeros([ntypes, ntypes, nr])
    fmax = -np.inf
    rhomax = -np.inf

    for i, n1 in enumerate(atnames):
        for j, n2 in enumerate(atnames):
            if j > i:
                continue
            elif i == j:
                rhor[n1] = prof(n1, r)
                rhomax = np.max(rhor[n1]) if rhomax < np.max(rhor[n1]) else rhomax
                z2r[i, j, :] = r * pair(n1, n2, r)
            else:
                z2r[i, j, :] = r * pair(n1, n2, r)
    z2r = np.where(z2r, z2r, z2r.transpose((1, 0, 2)))
    if rhomax < 2.0 * rhoemax:
        rhomax = 2.0 * rhoemax
    if rhomax < 100.0:
        rhomax = 100.0
    rho = np.linspace(0.0, rhomax, num=nrho, dtype=np.double)
    drho = rho[1] - rho[0]
    for i, n1 in enumerate(atnames):
        Fr[n1] = embed(n1, rho)

    write_file(atnames, outfilename, Fr, rhor, z2r, nrho, drho, nr, dr, rc)

    return


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        raise SystemExit("User interruption.")
