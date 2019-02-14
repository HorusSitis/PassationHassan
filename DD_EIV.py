#################################################################################################
## Etape IV : Prédictions. Choisir les paramètres du problème à résoudre par le modèle réduit. ##
#################################################################################################

### ------------ Reproduire pour des étapes I à IV. ------------ ###

tol=1e-10

xinf=0.0
yinf=0.0
xsup=1.0
ysup=1.0

#determiner le domaine fixe pour interpoler la solution

dimension=2

class PeriodicBoundary(SubDomain):
 # Left boundary is "target domain" G
 def inside(self, x, on_boundary):
  return on_boundary and not(near(x[0],xsup,tol) or near(x[1],ysup,tol))## merci à Arnold Douglas
 # Map right boundary (H) to left boundary (G)
 def map(self, x, y):
  for i in range(dimension):
   if near(x[i],1.0,tol):
    y[i]=0.0
   else:
    y[i]=x[i]

# maillage et fonctions tests du domaine fixe

mesh_fixe=Mesh("maillages_per/2D/maillage_fixe2D_am.xml")

V_fixe=VectorFunctionSpace(mesh_fixe, 'P', 3, constrained_domain=PeriodicBoundary())

# Performances

import time

### ------------ Etapes reproduites : dépendances directes de Main3D ------------ ###

nb_modes=N_mor

## mention="..." ## affectation effectuée en préambule, voir Main2D.py

if typ_msh=='gms':
 mesh_name="maillages_per/2D/maillage_trou2D"+mention+"_"+str(int(round(100*r_nouv,2)))+".xml"
 print(mesh_name)
 mesh_nouv=Mesh(mesh_name)

V_nouv=VectorFunctionSpace(mesh_nouv, "P", 3, constrained_domain=PeriodicBoundary())

# --------------------- SE1 : projection de la base POD sur le nouveau domaine --------------------- #

## On initialise le temps de calcul ##

start=time.time()

## Taille du maillage du domaine fixe ##

nb_noeuds_fixe=V_fixe.dim()

## Chargement de la base POD complète

phi_name='Phi'+'_dim'+str(Nsnap)+'_'+config+'_'+geo_p+'_'+"res"+str(res)+'_'+ordo+'_'+computer

print(phi_name)

with sh.open(repertoire_parent+phi_name) as phi_loa:
 Phi_prime_v = phi_loa["maliste"]

## Création de la base POD tronquée, sous forme vectorielle

Phi_mor=Phi_prime_v[:,range(0,nb_modes)]

## Extrpolation des fonctions de la base POD pour former le modèle réduit défini sur V_nouv

nb_noeuds_nouv=V_nouv.dim()
Phi_nouv_v=np.zeros((nb_noeuds_nouv,nb_modes))

list_pod_nouv=[]
phi_nouv=Function(V_nouv)
phi_fixe=Function(V_fixe)

for n in range(0,nb_modes):
 phi_fixe.vector().set_local(Phi_mor[:,n])
 # extrapolation du snapshot au domaine fixe
 phi_fixe.set_allow_extrapolation(True)
 phi_n_nouv=interpolate(phi_fixe,V_nouv)
 # affichage des modes extrapolés
 plot(phi_n_nouv)
 if fig_todo=='aff':
  plt.show()
 else:
  plt.savefig("Figures2D/phi_nouv_"+str(n+1)+"_"+config+'_'+geo_p+".png")
 plt.close()
 # on range le vecteur de POD interpolée dans la matrice Phi_nouv_v
 Phi_nouv_v[:,n]=phi_n_nouv.vector().get_local()

## On enregistre et imprime le temps d'éxécution de SE1

end=time.time()

print('se1 faite ',end-start,' secondes')
#sys.exit()#-------------------------------------
# --------------------- SE2 : résolution du modèle réduit --------------------- #

## On réinitialise le temps de calcul ##

start=time.time()

## On écrit les deux tenseurs qui comportent les coefficients de l'équation du modèle réduit : ceux-ci dépendent des vecteurs de la base POD projetée

#from PO23D import *

Coeff=calc_Ab_2D(V_nouv,mesh_nouv,Phi_nouv_v,r_nouv,cen_snap_ray,nb_modes)
A=Coeff[0]
b=Coeff[1]

## On résoud le modèle réduit

a_nouv=np.linalg.solve(A.T,-b)

## On enregistre et imprime le temps d'éxécution de SE2

end=time.time()

print('se2 faite ',end-start,' secondes')
# --------------------- SE3 : calcul du nouveau champ de vecteurs, affichage --------------------- #

## On réinitialise le temps de calcul ##

start=time.time()

## On initialise et affiche le champ chi_nouv

chi_nouv_v=np.dot(Phi_nouv_v,a_nouv)
chi_nouv=Function(V_nouv)
chi_nouv.vector().set_local(chi_nouv_v)

plot(chi_nouv)#, linewidth=0.55)
if fig_todo=='aff':
 plt.show()
else:
 plt.savefig("Figures2D/solROM_"+config+'_'+geo_p+str(int(round(100*r_nouv,2)))+".png")
plt.close()

## Exploitation du champ ainsi obtenu
r=r_nouv
rho=r_nouv

# Affichage des valeurs et erreurs de la solution périodique, quelle que soit la configuration
#err_per_ind_01(chi_n,cen,r,npas_err)

err_per_gr(cen_snap_ray,r_nouv,chi_nouv,npas_err,fig_todo)
# Tenseur de diffusion homogénéisé
## Intégrale de chi sur le domaine fluide
T_chi=np.zeros((2,2))
for k in range(0,2):
 for l in range(0,2):
  T_chi[k,l]=assemble(grad(chi_nouv)[k,l]*dx)
## Intégrale de l'identité sur le domaine fluide
por=(1-pi*r**2)
D=por*np.eye(2)
## Calcul et affichage du tenseur Dhom
Dhom_kMOR=D_k*(D+T_chi.T)
#print(('Tenseur Dhom_k',Dhom_k))
print('Coefficient Dhom_k11 '+conf_mess+', '+geo_mess+' valeur '+str(rho)+' MOR :',Dhom_kMOR[0,0])

## On enregistre et imprime le temps d'éxécution de SE3

end=time.time()

print('se3 faite ',end-start,' secondes')
#sys.exit()#-------------------------------------
# --------------------- SE4 : comparaison avec la méthode des éléments finis --------------------- #

## On réinitialise le temps de calcul ##

start=time.time()

## On réinitialise le champ chi_nouv pour la méthode des éléments finis

#res=20
chi_nouv=snapshot_circ_per(cen_snap_ray,r_nouv,res)

## Exploitation du champ ainsi obtenu
rho=r_nouv
r=r_nouv

# Affichage des valeurs et erreurs de la solution périodique, quelle que soit la configuration
#err_per_ind_01(chi_n,cen,r,npas_err)
err_per_gr(cen_snap_ray,r_nouv,chi_nouv,npas_err,fig_todo)

# Tenseur de diffusion homogénéisé
## Intégrale de chi sur le domaine fluide
T_chi=np.zeros((2,2))
for k in range(0,2):
 for l in range(0,2):
  T_chi[k,l]=assemble(grad(chi_nouv)[k,l]*dx)
## Intégrale de l'identité sur le domaine fluide
por=(1-pi*r**2)
D=por*np.eye(2)
print('Noeuds :',V_nouv.dim())
print('Porosité :',por)
## Calcul et affichage du tenseur Dhom
Dhom_kMEF=D_k*(D+T_chi.T)
#print(('Tenseur Dhom_k',Dhom_k))
print('Coefficient Dhom_k11 '+conf_mess+', '+geo_mess+' valeur '+str(rho)+ ' MEF :',Dhom_kMEF[0,0])

## Comparaison

err_rel=100*(Dhom_kMOR[0,0]-Dhom_kMEF[0,0])/Dhom_kMEF[0,0]
print('Erreur relative MEF-MOR :', err_rel , ' pourcent')

## Sortie graphique

plot(chi_nouv)#, linewidth=0.55)
if fig_todo=='aff':
 plt.show()
else:
 plt.savefig("Figures2D/solFEM_"+config+'_'+geo_p+str(int(round(100*r_nouv,2)))+".png")
plt.close()

## On enregistre et imprime le temps d'éxécution de SE4

end=time.time()

print('se4 faite ',end-start,' secondes')
