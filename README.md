# Kohn-Sham inversion code and QMC results

This repository contains the implementation of the Kohn-Sham inversion algorithm and the inversion results discussed in [Accurate Kohn-Sham auxiliary system from the ground state density of solids, arxiv:2207.03919, A. Aouina, M. Gatti, S. Chen, S. Zhang and L. Reining.](https://arxiv.org/abs/2207.03919)

# Results
The inverted exchange-correlation potentials obtained for the QMC densities are in the files ``si_qmc_dens_vxc.csv`` and ``nacl_qmc_dens_vxc.csv``. They are expressed in atomic units.  


# Silicon 
We considered the primitive cell, which is a Face Centered Cubic (FCC) cell: it contains two
Si atoms located at $(0, 0, 0)$ and $(\frac14, \frac14 , \frac14 )$. The volume of the primitive cell is
a quarter of the conventional cubic one which has a lattice constant of $10.263$
Bohr. To construct the KS potential, we used a Ne-core pseudopotential. We used a plane wave
cutoff of $25$ Ry, and $6 × 6 × 6$ k-point grid.

# NaCl 
We used the primitive
cell which contains one Na atom at lattice points and one Cl atom at bulk
center. The volume of this cell is one quarter the FCC cubic one. The lattice
constant used for the latter is $10.7563$ Bohr. The plane wave cutoff is fixed
to $40$ Ry and we used $6 × 6 × 6$ k-point grid. In this case a He-core pseudo-
potential is used for Na and Ne-core one for Cl. This implies eight occupied
valence bands. 

# Code
The inversion code is obtained by making a small modification to the [PyKSolver](https://github.com/aouinaayoub/PyKSolver) code which is used to solve the Kohn-Sham equations. 

To try the code for Silicon:
- Download the non-local part of the pseudo-potential of Silicon from here: [vnl.dat](https://zenodo.org/record/7661254/files/vnl.tar.gz?download=1) and put it in the same directory. 
- Then, execute the `pyksinvert.py` file:
```bash 
python pyksinvert.py
``` 