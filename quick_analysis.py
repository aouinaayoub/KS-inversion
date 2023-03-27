import numpy as np
import matplotlib.pyplot as plt 
from input_inv import * 
import sys 
####
n_rgrid=n_rgrid1 
atom2_pos=0
p1=np.round(np.dot(np.linalg.inv(MtR)*a_l,[atom2_pos,atom2_pos,atom2_pos])*n_rgrid)[0]
p2xy=np.round(np.dot(np.linalg.inv(MtR)*a_l,[atom2_pos,atom2_pos,1+atom2_pos])*n_rgrid)[0]
p2z=np.round(np.dot(np.linalg.inv(MtR)*a_l,[atom2_pos,atom2_pos,1+atom2_pos])*n_rgrid)[-1]
p1,p2xy,p2z = int(p1),int(p2xy),int(p2z)
def get_dens_inSD(dens): 
  dens3D=np.reshape(dens,(n_rgrid,n_rgrid,n_rgrid))
  diag=[]
  ### 001 --> 11-1
  for i in range(n_rgrid+1):
    ii=(p1+i)%n_rgrid
    jj=(p1+i)%n_rgrid
    kk=(p1-i)%n_rgrid
    diag.append(dens3D[ii,jj,kk])
    #ind_dens.append(ii+ jj*n_rgrid+ kk * n_rgrid**2)
  ### 110 --> 001  
  for i in range(2,n_rgrid*2,2):
    ii=p2xy %n_rgrid
    jj=p2xy %n_rgrid
    kk=(i+p2z)%n_rgrid
    diag.append(dens3D[ii,jj,kk])
  #111
  for i in range(p1,n_rgrid+1+p1):
    ii=i%n_rgrid
    jj=ii
    kk=ii
    diag.append(dens3D[ii,jj,kk])
  return np.array(diag)

def analyze_ref_inversion(file_name, get_variables=0, max_it=-1,save=0, save_name="quick_analysis"): 
    """ The function returns: error_dens_atoms, error_max_it_pot_atom, error_max_it, vxc_atoms_it,error_max_it_pot, vxc_it[max_it]"""
    
    loaded_file=np.load(file_name)
    dens_it, vxc_it= loaded_file["dens_history"].real ,loaded_file["potxc_history"].real
    dens_it,vxc_it=dens_it[:max_it], vxc_it[:max_it]
    try : 
        dens_ref= np.genfromtxt(sys.argv[2]) 
    except: 
        dens_ref= dens_it[max_it]
    
    vxc_ref=vxc_it[max_it]
    error_max_it=[np.max(dens_it[i]/dens_ref - 1)*100 for i in range(len(dens_it))]
    error_mean_it=[np.mean(dens_it[i]/dens_ref - 1)*100 for i in range(len(dens_it))]
    
    shift_vxc_it=[np.mean(vxc_ref)-np.mean(vxc_it[i]) for i in range(len(vxc_it))]

    print(" > density max error of the last it is {:.5f} % " .format(error_max_it[-1]))
    
    print(" > the minimum of density max errors is {:.5f} %  at iteration {}" .format(np.min(error_max_it[1:]), 1+np.argmin(error_max_it[1:])))

    f,sub=plt.subplots(2,3,figsize=(20,10))    
    sub[0,0].plot(error_max_it[1:], label="error max")
    sub[0,0].plot(error_mean_it[1:], label="error mean")
    sub[0,0].legend()
    sub[0,0].set_xlabel("iteration")
    sub[0,0].set_ylabel("error on density (%)")
   
    sub[0,1].plot(get_dens_inSD(dens_it[-1]/dens_ref -1) *100)
    sub[0,1].set_ylabel("error on density (%)")
    sub[0,1].set_xlabel("points along the route")
   
    sub[0,2].plot(get_dens_inSD( (vxc_it[-1]+shift_vxc_it[-1])/vxc_ref -1) *100, label="iteration -1 (last)")
    sub[0,2].plot(get_dens_inSD( (vxc_it[-2]+shift_vxc_it[-2])/vxc_ref -1) *100, ls="dashdot", label="iteration -2 ")
    sub[0,2].plot(get_dens_inSD( (vxc_it[-3]+shift_vxc_it[-3])/vxc_ref -1) *100, ls="dashed", label="iteration -3")
    sub[0,2].set_ylabel("error on vxc (%)")
    sub[0,2].legend()
    sub[0,2].set_xlabel("points along the route")
    #plt.show()
   
    sub[1,1].plot(get_dens_inSD( vxc_it[-1]+shift_vxc_it[-1]), label="iteration -1 (last)")
    sub[1,1].plot(get_dens_inSD( vxc_it[-2]+shift_vxc_it[-2]), ls="dashdot", label="iteration -2 ")
    sub[1,1].plot(get_dens_inSD( vxc_it[-3]+shift_vxc_it[-3]), ls="dashed", label="iteration -3")
    sub[1,1].scatter(range(len(get_dens_inSD(vxc_ref))),get_dens_inSD( vxc_ref), label="ref", color="grey")
    sub[1,1].set_ylabel("vxc")
    sub[1,1].legend()
    sub[1,1].set_xlabel("points along the route")


    # max error of the potential vs error of the potential on atoms
    error_max_it_pot= [np.max(np.abs((vxc_it[i]+shift_vxc_it[i])/vxc_ref - 1))*100 for i in range(len(vxc_it))]
    sub[1,0].plot(error_max_it_pot, label="error max")
    sub[1,0].set_xlabel("iterations")
    sub[1,0].set_ylabel("error on potential (%)")
    sub[1,0].legend()
    #plt.show() 

    sub[1,2].scatter(dens_ref,vxc_it[-1]+shift_vxc_it[-1], label="iteration -1 (last one)")
    sub[1,2].scatter(dens_ref,vxc_it[-2]+shift_vxc_it[-2], label="iteration -2", marker="^")
    sub[1,2].scatter(dens_ref,vxc_it[-3]+shift_vxc_it[-3], label="iteration -3", marker="+")
    sub[1,2].set_xlabel("n(r)")
    sub[1,2].set_ylabel("vxc(r)")
    sub[1,2].legend()
    if save: 
        plt.savefig(save_name+".pdf")
    plt.show() 

    if get_variables:
        return  error_max_it, error_max_it_pot, vxc_it[max_it]

analyze_ref_inversion(sys.argv[1],save=1)
