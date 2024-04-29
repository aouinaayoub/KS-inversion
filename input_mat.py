import os
import numpy as np 

#-------------------------------------------------------------------------------------------------------------------------#
# Reference  density
densR_ref= np.genfromtxt("dens_ref_lda.dat") 

densR= np.copy(densR_ref)

# choice of the xc functional 
xc_type="LDA"

# choose output file name   
output_filename=xc_type 
# continue from last run 
continue_from_last_run=0 
last_run="invlda.npz" 

#-------------------------------------------------------------------------------------------------------------------------#

# parameters of the material in bohr 
a_l=np.array([10.263087]*3) 

# Volume of unit cell 
V_unitcell= np.dot(np.cross([a_l[0],0,0],[0,a_l[1],0]) , [0, 0,a_l[2]])/4 

# change of basis from conventional to primitive cell 
MtR=np.array([[0,a_l[0]/2,a_l[0]/2],[a_l[0]/2, 0,a_l[0]/2],[a_l[0]/2, a_l[0]/2,0]])

# number of occupied bands for semi-conductor and insulator
n_occup=4

# Grid of k 
kx=np.arange(-2,4)*1/6  #  sampling of K points 

# Cut-off energy 
Ecut=12.5

# nonlocal and local potenial directory 
origin_dir= "./"

# wrap up
input_dict={"a_l":a_l, "V_unitcell":V_unitcell, "MtR":MtR, "n_occup":n_occup, "kx":kx, "Ecut":Ecut, 
            "path_VextR": origin_dir+"VPS_loc.dat", "path_vnl": origin_dir+'vnl.dat', 
            "densR": densR, "densR_ref": densR_ref , "output": output_filename, "xc_type":xc_type, "dict_vnl": True, "continue_from_last_run":continue_from_last_run, "last_run":last_run} 
