#####################################################################################################################################
######################################### Etape II : extrapolation des clichés, domaine_fixe ########################################
#####################################################################################################################################


## chargement du snapshot pour l'indice courant
# Extrapolation au domaine Omega_fixe : aucune inclusion, khi défini sur [0,1]times[0,1]
chi.set_allow_extrapolation(True)
chi_fixe=interpolate(khi,V_fixe)##rapide




