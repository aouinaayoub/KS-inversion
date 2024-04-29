# Kohn-Sham Inversion Code and QMC Results

This repository contains the implementation of the Kohn-Sham inversion algorithm and the inversion results discussed in the paper "[Accurate Kohn-Sham auxiliary system from the ground state density of solids](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.107.195123)" by A. Aouina, M. Gatti, S. Chen, S. Zhang, and L. Reining.

## Results

The inverted exchange-correlation potentials obtained from the Quantum Monte Carlo (QMC) densities are provided in the following files:

- `si_qmc_dens_vxc.csv`
- `nacl_qmc_dens_vxc.csv`

The potentials are expressed in atomic units.

## Silicon

For silicon, we considered the primitive cell, which is a Face Centered Cubic (FCC) cell containing two Si atoms located at coordinates (0, 0, 0) and (1/4, 1/4, 1/4). The volume of the primitive cell is one quarter of the conventional cubic one, with a lattice constant of 10.263087 Bohr. In constructing the Kohn-Sham potential, we used a [Ne-core pseudopotential](https://github.com/aouinaayoub/KS-inversion/blob/main/pseudo-potentials/si.upf) and a plane wave cutoff of 25 Ry. Additionally, a 6 × 6 × 6 k-point grid was employed.

## NaCl

For NaCl, we used the primitive cell containing one Na atom at lattice points and one Cl atom at the bulk center. The volume of this cell is one quarter of the FCC cubic one. The lattice constant used is 10.7563 Bohr. The plane wave cutoff was set to 40 Ry, and a 6 × 6 × 6 k-point grid was used. In this case, a [He-core pseudopotential](https://github.com/aouinaayoub/KS-inversion/blob/main/pseudo-potentials/Na.upf) was employed for Na and a [Ne-core pseudopotential](https://github.com/aouinaayoub/KS-inversion/blob/main/pseudo-potentials/Cl.upf) for Cl, resulting in eight occupied valence bands.

## Code

The inversion algorithm is implemented as an *xc functional class*  in the [PyKSolver](https://github.com/aouinaayoub/PyKSolver) code, which is used to solve the Kohn-Sham equations. 

To run the code for silicon:

1. Download the non-local part of the pseudo-potential of silicon from [here](https://zenodo.org/record/7661254/files/vnl.tar.gz?download=1) and place it in the same directory.
2. Execute the `pyksolver.py` file:
   
   ```bash 
   python pyksolver.py
   ```

Note: Make sure you have the necessary dependencies and libraries installed before running the code.

Please refer to the original paper for more details on the Kohn-Sham inversion algorithm and the obtained results. 
