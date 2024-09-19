from pyscf import gto, scf
import sys
import os
import numpy as np
import fortranformat as ff

basis = 'ccpvdz'
ifile = sys.argv[2]
ofile = sys.argv[3]

def read_gauss_input(ifile):
    with open(ifile, "r") as ifile:
        gau_input = ifile.readlines()

    tokens = gau_input[0].split()

    natoms = int(tokens[0])
    chrg = int(tokens[2])
    spin = int(tokens[3]) -1 
    coords = np.empty((natoms, 4))
    atomtypes=[]
    for i, line in enumerate(gau_input[1 : 1 + natoms]):
        line = line.split()
        line_re = line[:-1]
        coords[i] = line_re
    return natoms, chrg, spin, coords

def convert_to_xyz(data):
    coord_empty = ''
    coord_final =''
    for i, value in enumerate(data):
        coord_empty = ' '.join(f'{x:.8f}' if not float(x).is_integer() else f'{int(x)}' for x in value) 
        coord_final += coord_empty
        if i < len(data)-1:
            coord_final += ';'
    return coord_final

def pyscf_inp():
    mol=gto.Mole()
    mol.spin = s
    mol.charge = c
    mol.basis = basis
    mol.atom = cr
    mol.unit = 'bohr'
    mol.verbose = 4
    mol.build()
    mf = scf.RHF(mol)
    scf_ener = mf.scf()
    nuc_ener = mf.energy_nuc()
    g = mf.nuc_grad_method()
    grad = g.kernel()
    return grad, scf_ener, nuc_ener


def read_gauss_out(ofile, energy, natoms, gradient):

    headformat = ff.FortranRecordWriter("4D20.12")
    bodyformat = ff.FortranRecordWriter("3D20.12")

    f = open(ofile, "w")
    head = [energy, 0 , 0, 0]
    headstring = headformat.write(head)
    f.write(headstring + "\n")

    for i in range(natoms):
        output = bodyformat.write(gradient[i])
        f.write(output + "\n")

    # polarizability and dipole derivatives are set to zero
    polarizability = np.zeros((2, 3))
    dipole_derivative = np.zeros((3 * natoms, 3))

    for i in range(2):
        output = bodyformat.write(polarizability[i])
        f.write(output + "\n")

    for i in range(3 * natoms):
        output = bodyformat.write(dipole_derivative[i])
        f.write(output + "\n")

    f.close()


if __name__=="__main__":
    n , c, s, coord = read_gauss_input(ifile)
    cr = convert_to_xyz(coord)
    gradient, energy_scf, energy_nuc = pyscf_inp()
    read_gauss_out(ofile, energy_scf, n, gradient)
