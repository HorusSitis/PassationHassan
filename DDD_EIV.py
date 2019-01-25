#################################################################################################
## Etape IV : Prédictions. Choisir les paramètres du problème à résoudre par le modèle réduit. ##
#################################################################################################

### ------------ Reproduire éventuellement pour des étapes ultérieures. Laisser seulement dans DDD_fun_obj ? ------------ ###

tol=1e-10

xinf=0.0
yinf=0.0
zinf=0.0
xsup=1.0
ysup=1.0
zsup=1.0

dimension=3

class PeriodicBoundary(SubDomain):
 # Left boundary is "target domain" G
 def inside(self, x, on_boundary):
  return on_boundary and not(near(x[0],xsup,tol) or near(x[1],ysup,tol) or near(x[2],zsup,tol))
 # Map right boundary (H) to left boundary (G)
 def map(self, x, y):
  for i in range(dimension):
   if near(x[i],1.0,tol):
    y[i]=0.0
   else:
    y[i]=x[i]

# maillage du domaine fixe

if dom_fixe=='':
 mesh_fixe=Mesh("maillages_per/3D/cubesphere_periodique_triangle.xml")
elif dom_fixe=='0001':
 mesh_fixe=Mesh("maillages_per/3D/cubesphere_periodique_triangle_0001fixe.xml")
elif dom_fixe=='0000':
 mesh_fixe=Mesh("maillages_per/3D/cubesphere_periodique_triangle_0000fixe.xml")

# fonctions test du domaine fixe

V_fixe=VectorFunctionSpace(mesh_fixe,'P',2,constrained_domain=PeriodicBoundary())

### ------------ Etapes reproduites : dépendances directes de Main3D ------------ ###

#r_nouv=0.22
nb_modes=N_mor

if typ_msh=='gms':
 mesh_nouv=Mesh("maillages_per/3D/cubesphere_periodique_triangle_"+str(int(round(100*r_nouv,2)))+".xml")
#else:mesh_nouv=creer_maill_sph(cen,r_nouv,res)

V_nouv=VectorFunctionSpace(mesh_nouv, "P", 2, constrained_domain=PeriodicBoundary())

# --------------------- SE1 : projection de la base POD sur le nouveau domaine --------------------- #

nb_noeuds_fixe=V_fixe.dim()

## Chargement de la base POD complète

phi_name='Phi'+dom_fixe+'_dim'+str(Nsnap)+'_'+config+'_'+geo_p+'_'+ordo+'_'+computer

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
 # on range le vecteur de POD interpolée dans la matrice Phi_nouv_v
 Phi_nouv_v[:,n]=phi_n_nouv.vector().get_local()

## Stockage de la matrice du modèle réduit



print('se1 faite')
#sys.exit()#-------------------------------------
# --------------------- SE2 : résolution du modèle réduit --------------------- #

## On écrit les deux tenseurs qui comportent les coefficients de l'équation du modèle réduit : ceux-ci dépendent des vecteurs de la base POD projetée

#from PO23D import *

Coeff=calc_Ab(V_nouv,mesh_nouv,Phi_nouv_v,r_nouv,cen_snap_ray,nb_modes)
A=Coeff[0]
b=Coeff[1]

## On résoud le modèle réduit

a_nouv=np.linalg.solve(A.T,-b)

print('se2 faite')
# --------------------- SE3 : calcul du nouveau champ de vecteurs, affichage --------------------- #

chi_nouv_v=np.dot(Phi_nouv_v,a_nouv)
chi_nouv=Function(V_nouv)
chi_nouv.vector().set_local(chi_nouv_v)

plot(chi_nouv, linewidth=0.55)
if fig_todo=='aff':
 plt.show()
#else:
plt.close()

## Exploitation du champ ainsi obtenu
r=r_nouv
rho=r_nouv

# Affichage des valeurs et erreurs de la solution périodique, quelle que soit la configuration
#err_per_ind_01(chi_n,cen,r,npas_err)

err_per_gr(cen_snap_ray,r_nouv,chi_nouv,npas_err,fig_todo)
# Tenseur de diffusion homogénéisé
## Intégrale de chi sur le domaine fluide
T_chi=np.zeros((3,3))
for k in range(0,3):
 for l in range(0,3):
  T_chi[k,l]=assemble(grad(chi_nouv)[k,l]*dx)
## Intégrale de l'identité sur le domaine fluide
if config=='sph_un':
 D=(1-4/3*pi*r**3)*np.eye(3)
elif config=='cyl_un':
 D=(1-pi*r**2)*np.eye(3)
else :
 D=(1-4/3*pi*r_s**3-pi*r_c**2)*np.eye(3)
## Calcul et affichage du tenseur Dhom
Dhom_kMOR=D_k*(D+T_chi.T)
#print(('Tenseur Dhom_k',Dhom_k))
print('Coefficient Dhom_k11 '+conf_mess+', '+geo_mess+' variable valeur '+str(rho)+' MOR :',Dhom_kMOR[0,0])

print('se3 faite')
#sys.exit()#-------------------------------------
# --------------------- SE4 : comparaison avec la méthode des éléments finis --------------------- #

res=20
if config=='sph_un':
 chi_nouv=snapshot_sph_per(cen_snap_ray,r_nouv,res)
elif config=='cyl_un':
 chi_nouv=snapshot_cyl_per(cen_snap_ray,r_nouv,res)


## Exploitation du champ ainsi obtenu
rho=r_nouv
r=r_nouv

# Affichage des valeurs et erreurs de la solution périodique, quelle que soit la configuration
#err_per_ind_01(chi_n,cen,r,npas_err)
err_per_gr(cen_snap_ray,r_nouv,chi_nouv,npas_err,fig_todo)

# Tenseur de diffusion homogénéisé
## Intégrale de chi sur le domaine fluide
T_chi=np.zeros((3,3))
for k in range(0,3):
 for l in range(0,3):
  T_chi[k,l]=assemble(grad(chi_nouv)[k,l]*dx)
## Intégrale de l'identité sur le domaine fluide
if config=='sph_un':#'sphère unique':
 por=(1-4/3*pi*r**3)
 D=por*np.eye(3)
elif config=='cyl_un':#'cylindre unique':
 por=(1-pi*r**2)
 D=por*np.eye(3)
else :
 por=(1-4/3*pi*r_s**3-pi*r_c**2)
 D=por*np.eye(3)
print('Noeuds :',V_nouv.dim())
print('Porosité :',por)
## Calcul et affichage du tenseur Dhom
Dhom_kMEF=D_k*(D+T_chi.T)
#print(('Tenseur Dhom_k',Dhom_k))
print('Coefficient Dhom_k11 '+conf_mess+', '+geo_mess+' variable valeur '+str(rho)+ ' MEF :',Dhom_kMEF[0,0])

## Comparaison

err_rel=100*(Dhom_kMOR[0,0]-Dhom_kMEF[0,0])/Dhom_kMEF[0,0]
print('Erreur relative MEF-MOR :', err_rel , ' pourcent')




