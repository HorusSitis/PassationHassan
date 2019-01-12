#####################################################################################################################################
######################################### Etape II : extrapolation des clichés, domaine_fixe ########################################
#####################################################################################################################################

#for i in range(1,1+Npas):
# r=i*0.2#05
# u=snapshot_circ_per([c_x,c_y],r,res)
# ## chargement du snapshot pour l'indice courant
# # Extrapolation au domaine Omega_fixe : aucune inclusion, khi défini sur [0,1]times[0,1]
# u.set_allow_extrapolation(True)
# u_fixe=interpolate(u,V_fixe)##rapide
# #u_fixe = project(u, V_fixe)##lent
# plot(u_fixe)
# plt.show()
# plt.close()

## On recharge les clichés stockés à l'étape I, aec le format hdf5
#if [c_x,c_y]==[0.5,0.5]:
# suffixe="inc_centre/"
#elif [c_x,c_y]==[0.0,0.0]:
# suffixe="coins/"

#repertoire_lecture=repertoire_parent+suffixe

#for i in range(1,1+Npas):
# r=rayons[i]
# # Création de l'espace dans lequel vit le cliché à charger
# mesh_c_r=creer_maill_circ([c_x,c_y],r,res)
# V=VectorFunctionSpace(mesh_c_r, 'P', 2, constrained_domain=PeriodicBoundary())
# chi_i=Function(V)
# # Chargement du cliché, réalisé à l'étape I
# lecture_fichiers_restart_hdf5(chi_i,repertoire_lecture,r,mesh_c_r)








# Matrice des snapshots : dans l'espace des fonctions extrapolées


