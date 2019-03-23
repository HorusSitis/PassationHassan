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

if dom_fixe=="am":
 mesh_fixe=Mesh("maillages_per/2D/maillage_fixe2D_am.xml")
elif dom_fixe=="multiray":
 mesh_fixe=Mesh("maillages_per/2D/maillage_fixe2d_"+dom_fixe+".xml")
elif config=='compl':
 mesh_fixe=Mesh("maillages_per/2D/maillage_trous2D_"+geo_p+"_fixe.xml")
elif dom_fixe=="ray_min":
 if config=='cer_un':
  mesh_fixe=Mesh('maillages_per/2D/maillage_trou2D_5.xml')

V_fixe=VectorFunctionSpace(mesh_fixe, 'P', VFS_degree, constrained_domain=PeriodicBoundary())

# Performances

import time

## mention="..." ## affectation effectuée en préambule, voir Main2D.py

### ------------ Etapes spécifiques à la construction du modèle réduit ------------ ###

nb_modes=N_mor

r=r_nouv
# Création du domaine d'intégration sur le maillage fixe : pour le calcul de Dhom
class DomPhysFluide(SubDomain):
 def inside(self, x, on_boundary):
  return True if ((x[0]-0.5)**2+(x[1]-0.5)**2>=r**2) else False
# marquage du domaine
dom_courant=DomPhysFluide()
subdomains=MeshFunction('size_t',mesh_fixe,mesh_fixe.topology().dim())
subdomains.set_all(1)
dom_courant.mark(subdomains,12829)
dxf=Measure("dx", domain=mesh_fixe, subdomain_data=subdomains)

if typ_msh=='gms':
 mesh_repository="maillages_per/2D/"
 if config!='compl':
  mesh_name="maillage_trou2D"+mention+"_"+str(int(round(100*r_nouv,2)))+".xml"
 else:
  mesh_name="maillage_trous2D_"+geo_p+"_"+str(int(round(100*r_nouv,2)))+".xml"

print("Maillage de Omega_nouv",mesh_repository+mesh_name)
mesh_nouv=Mesh(mesh_repository+mesh_name)

V_nouv=VectorFunctionSpace(mesh_nouv, "P", VFS_degree, constrained_domain=PeriodicBoundary())

## Taille du maillage du domaine fixe ##

nb_noeuds_fixe=V_fixe.dim()

## Chargement de la base POD complète

phi_name='Phi'+dom_fixe+'_dim'+str(Nsnap)+'_'+config+'_'+geo_p+'_deg'+str(VFS_degree)+'_'+"res"+str(res)+'_'+ordo+'_'+computer

print(phi_name)

with sh.open(repertoire_parent+phi_name) as phi_loa:
 Phi_prime_v = phi_loa["maliste"]

## Création de la base POD tronquée, sous forme vectorielle

Phi_mor=Phi_prime_v[:,range(0,nb_modes)]

### ------------ Fin ------------ ###

# --------------------- SE1 : projection de géométrie courante sur le domaine fixe, création de maillage --------------------- #

## On initialise le temps de calcul ##

start_se1=time.time()

# Initialisation des fonctions POD

phi_fixe=Function(V_fixe)

# Raffinement du maillage : pour le moment, cas d'une inclusion unique centrée
start_refi=time.time()
mesh_r_fixe=mesh_fixe
r=r_nouv
for i in range(1,1+Nrefine):
 print('raffinement',i)
 markers = MeshFunction("bool", mesh_fixe, mesh_fixe.topology().dim())
 markers.set_all(False)
 for c in cells(mesh_fixe):
  for f in facets(c):
   if ((f.midpoint()[0]-cen_snap_ray[0])**2+(f.midpoint()[1]-cen_snap_ray[1])**2<=(r*(1+crow/i))**2) and ((f.midpoint()[0]-cen_snap_ray[0])**2+(f.midpoint()[1]-cen_snap_ray[1])**2>=(r*(1-crow/i))**2):
    markers[c]=True
  for v in vertices(c):
   if (v.point().x()-cen_snap_ray[0])**2+(v.point().y()-cen_snap_ray[1])**2>=(r*(1-crow/i))**2 and (v.point().x()-cen_snap_ray[0])**2+(v.point().y()-cen_snap_ray[1])**2<=(r*(1+crow/i))**2:
    markers[c]=True
 mesh_r_fixe=refine(mesh_r_fixe, markers, redistribute=True)
end=time.time()
#print('Raffinemenent du maillage :',end-start,'secondes')
tps_refi=end-start_refi

if fig_todo=='aff':
 plot(mesh_r_fixe)
 plt.title("Maillage raffiné autour de la frontière physique, rayon "+str(int(round(100*r_nouv,2))))
 plt.show()
 plt.close()

V_r_fixe=VectorFunctionSpace(mesh_r_fixe, "P", VFS_degree, constrained_domain=PeriodicBoundary())

nb_noeuds_r_fixe=V_r_fixe.dim()
Phi_r_fixe_v=np.zeros((nb_noeuds_r_fixe,nb_modes))

# Interpolation sur le maillage raffiné
start_interp=time.time()
if Nrefine>0:
 for n in range(0,nb_modes):
  phi_fixe.vector().set_local(Phi_mor[:,n])
  # extrapolation du snapshot au domaine fixe
  phi_fixe.set_allow_extrapolation(True)
  phi_r_fixe=interpolate(phi_fixe,V_r_fixe)
  # affichage des modes extrapolés
  plot(phi_r_fixe)
  plt.title("Phi "+str(n+1)+" sur le maillage raffiné",fontsize=30)
  if fig_todo=='aff':
   plt.show()
  elif fig_todo=='save':
   plt.savefig("Figures2D/phi_fixe_r_"+str(n+1)+"_"+config+'_'+geo_p+".png")
  plt.close()
  # on range le vecteur de POD interpolée dans la matrice Phi_nouv_v
  Phi_r_fixe_v[:,n]=phi_r_fixe.vector().get_local()
 # on mesure le temps d'interpolation
 end=time.time()
 tps_interp=end-start_interp

## On enregistre et imprime le temps d'éxécution de SE1

end=time.time()

print('raffinement fait',tps_refi,'secondes')
print('interpolation faite',tps_interp,'secondes')
print('se1 faite',end-start_se1,'secondes')
sys.exit()#-------------------------------------
# --------------------- SE2 : chargement de a base POD et résolution du modèle réduit --------------------- #

## On réinitialise le temps de calcul ##

start=time.time()

## -------------------------------------------------------------------------------------- ##
## Extrapolation des fonctions de la base POD pour former le modèle réduit défini sur V_nouv

#nb_noeuds_nouv=V_nouv.dim()
#Phi_nouv_v=np.zeros((nb_noeuds_nouv,nb_modes))

#list_pod_nouv=[]
#phi_nouv=Function(V_nouv)
phi_fixe=Function(V_fixe)

## -------------------------------------------------------------------------------------- ##
















## On écrit les deux tenseurs qui comportent les coefficients de l'équation du modèle réduit : ceux-ci dépendent des vecteurs de la base POD fixe resreinte au domaine fluide virtuel

if test_snap=='i_per':
 Coeff=calc_Ab_2D(V_nouv,mesh_nouv,Phi_nouv_v,r_nouv,cen_snap_ray,nb_modes)
else:
 Coeff=calc_Ab_compl(V_nouv,mesh_nouv,Phi_nouv_v,nb_modes,test_snap)
 print('modèle réduit complexe utilisé')

A=Coeff[0]
b=Coeff[1]

## On résoud le modèle réduit

a_nouv=np.linalg.solve(A.T,-b)

## On enregistre et imprime le temps d'éxécution de SE2

end=time.time()

print('se2 faite ',end-start,' secondes')
# --------------------- SE3 : calcul du tenseur homogénéisé sur le domaine fixe, affichage de chi_nouv_fixe restreint au domaine fluide virtuel --------------------- #

## On réinitialise le temps de calcul ##

start=time.time()

## On initialise et affiche le champ chi_nouv

chi_nouv_v=np.dot(Phi_nouv_v,a_nouv)
chi_nouv=Function(V_nouv)
chi_nouv.vector().set_local(chi_nouv_v)

plot(chi_nouv)#, linewidth=0.55)
plt.title("Solution ROM", fontsize=30)#"Rho = 0,"+str(int(round(100*r_nouv,2))),fontsize=40)
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

#err_per_gr(cen_snap_ray,r_nouv,chi_nouv,npas_err,fig_todo)
# Tenseur de diffusion homogénéisé
## Intégrale de chi sur le domaine fluide
T_chi=np.zeros((2,2))
for k in range(0,2):
 for l in range(0,2):
  T_chi[k,l]=assemble(grad(chi_nouv)[k,l]*dx)
## Intégrale de l'identité sur le domaine fluide
if config!='compl':
 por=(1-pi*r**2)
else:
 por=1-pi*(r**2+0.15**2)
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
if test_snap=='i_per':
 chi_nouv=snapshot_circ_per(cen_snap_ray,r_nouv,res)
else:
 chi_nouv=snapshot_compl_per(geo_p,r_nouv,cen_snap_ray,mention,test_snap)

## Exploitation du champ ainsi obtenu
rho=r_nouv
r=r_nouv

# Affichage des valeurs et erreurs de la solution périodique, quelle que soit la configuration
#err_per_ind_01(chi_n,cen,r,npas_err)
for npas_test in []:#30,40]:#7,15,16,60,125,250]:
 err_per_gr(cen_snap_ray,r_nouv,chi_nouv,npas_test,fig_todo)

# Tenseur de diffusion homogénéisé
## Intégrale de chi sur le domaine fluide
T_chi=np.zeros((2,2))
for k in range(0,2):
 for l in range(0,2):
  T_chi[k,l]=assemble(grad(chi_nouv)[k,l]*dx)
## Intégrale de l'identité sur le domaine fluide
if config!='compl':
 por=(1-pi*r**2)
else:
 por=1-pi*(r**2+0.15**2)
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
plt.title("Solution EF",fontsize=30)
if fig_todo=='aff':
 plt.show()
else:
 plt.savefig("Figures2D/solFEM_"+config+'_'+geo_p+str(int(round(100*r_nouv,2)))+".png")
plt.close()

## On enregistre et imprime le temps d'éxécution de SE4

end=time.time()

print('se4 faite ',end-start,' secondes')
