import numpy as np
import os

#-------------------------------------------------------------------------------------------------------------------------#

# choose output file name   
output_filename="invlda"

# The reference density
densR_ref= np.genfromtxt("dens_ref_lda.dat") 

# intial guess for the xc pot  
potxcR= 0.3* (-(3./np.pi)**(1./3.)*densR_ref**(1./3.)) -0.2 # start with 0.3 vlda_x 

# continue from last run 
continue_from_last_run=1 
last_run="invlda.npz" 

#-------------------------------------------------------------------------------------------------------------------------#

# parameters of the material in bohr 
a_l=np.array([10.263087]*3) 

# change of basis from conventional to primitive cell 
MtR=(np.array([[0,a_l[0]/2,a_l[0]/2],[a_l[0]/2, 0,a_l[0]/2],[a_l[0]/2, a_l[0]/2,0]]))

#number of occupied bands for semi-conductor and insulator
n_occup=4

#Volume of the primitive cell ( for Si is 1/4 the volume of the conventional cell)  
V_unitcell= np.dot(np.cross([a_l[0],0,0],[0,a_l[1],0]), [0, 0,a_l[2]])/4  

# Grid of k 
kx=np.arange(-2,4)*1/6  #  sampling of K points 

# Cut-off energy 
Ecut=12.5

# local potential  
V_extR =np.genfromtxt("VPS_loc.dat")

# grid in real space, it depends on the cutoff 
n_rgrid1,n_rgrid2,n_rgrid3 = [round(len(V_extR)**(1/3))]*3

# nonlocal potenial 
origin_dir= "."
path_vnl=origin_dir+'vnl.dat'

# construct or load matrix 
matrix_dir="./matrix"
flagmatrix=0 
if not os.path.exists(matrix_dir): 
  os.makedirs(matrix_dir)
  flagmatrix=1
## create a matrix for each n_pkgmat 
n_pkgmat=1

# damping parameter 
alpha=5e-5
