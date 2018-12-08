import sys
from dolfin import *
import pylab as py

########################################################################
### Commandes de Cyrille pour stocker et charger un maillage ###########
#ECRITURE                                                              #
#                                                                      #
#Name_mesh="%s/mesh_cylindre_fixe" %(repertoire_ecriture)+".xdmf"      #
#File(Name_mesh) << mesh                                               #
#                                                                      #
#LECTURE                                                               #
#                                                                      #
#Name_mesh="%s/mesh_cylindre_t=" %(repertoire_init)+str(T_init)+".xdmf"#
#mesh=Mesh(Name_mesh)                                                  #
########################################################################

def creation_fichier_pourecriture_champ_hdf5(repertoire_final,mesh):
 "Creation des dichiers pour ecrire la solution"
 kh_file = File("%s/primi_hom.pvd" %(repertoire_final))                       
 KH_SAVE=HDF5File(mesh.mpi_comm(),"%s/khi_save.hdf5" %(repertoire_final), "w")
 return kh_file,KH_SAVE

def ecriture_champ_hdf5(kh_file,KH_SAVE,sol_n,kfic,file_rayon_ecriture,r,cen,res):
 c_x,c_y=cen[0],cen[1]
 "Ecriture des differents champs dans des fichiers"
 kh_file << sol_n
 KH_SAVE.write(sol_n,"khi",kfic)                                           
 #file_rayon_ecriture.write(str(kfic)+"\t"+str(r)+"\n"+str(res))
 return
