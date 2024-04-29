
from input_mat import input_dict        
from utils import *
import pickle 

xc_type="inversion"
input_dict["output"]="si_" + xc_type 

si=system(input_dict)
     
ham_instance = Hamiltonian_template(si)
# with open("si_ham_template","rb") as f :
#     #pickle.dump(ham_instance,f) 
#     ham_instance= pickle.load(f)

xc_inv = vxc_inverter()
xc_inv.xc_type =xc_type 
# initial guess for the potential 
xc_inv.potxcR = (-(3./np.pi)**(1./3.)*si.densR_ref**(1./3.)) - .4

###
ks_solv = KS_solver(si) 

for n_it in range(1,999):
    xc_inv.update_xc(si) 
    ks_solv.get_dens_parallel(si,xc_inv,ham_instance, mixing=0.75, njobs=1)
    reporting(n_it,si, ks_solv, xc_inv)
