#######################################################################################################################################
## Etape I : réalisation des clichés, avec la méthode des éléments finis. Calcul du tenseur d'homogénéisation. Stockage dans snap2D/ ##
#######################################################################################################################################

### ------------ Reproduire éventuellement pour des étapes ultérieures. Laisser seulement dans DDD_fun_obj ? ------------ ###

tol=1e-10

xinf=0.0
yinf=0.0
zinf=0.0
xsup=1.0
ysup=1.0
zsup=1.0

#determiner le domaine fixe pour interpoler la solution

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

### Maillage sur le domaine Omega_fixe : aucune inclusion, porosité égale à 1 ###

domaine_fixe=Box(Point(xinf,yinf,zinf),Point(xsup,ysup,zsup))

if typ_msh=='gms':
 mesh_fixe=Mesh("maillages_per/3D/cubesphere_periodique_triangle_0001fixe.xml")
else:
 mesh_fixe=generate_mesh(domaine_fixe,res_fixe)

V_fixe=VectorFunctionSpace(mesh_fixe, "P", 2, constrained_domain=PeriodicBoundary())

#if fig_todo=='aff':
# #représentation graphique du maillage
# plot(mesh_fixe)
# plt.show()
# plt.close()
#else:
# #sauvegarde de la figure
# plot(mesh_fixe)
# plt.savefig("Figures3D/mesh_fixe.png")
# plt.close()




#chi_s=snapshot_sph_per([0.5,0.5,0.5],0.25,res)
#chi_c=snapshot_cyl_per([0.5,0.5],0.1,slices_cyl,res)
#plot(chi_s)
#if fig_todo=='aff':
# plt.show()
#plt.close()
#plot(chi_c)
#plt.show()
#plt.close()

##Snapshot unique, avec les maillages de Cyrille





## Boucle pour la création des snapshots, avec un paramètre pouvant être le rayon d'une inclusion circulaire, ou l'emplacement de son centre 

# Pour avoir des fonctions "top-level" à paralléliser

## Sphère unique

#if geo_p=='rayon':
cen_snap_ray=[0.5,0.5,0.5]
def snap_ray(r_par):
 chi_r=snapshot_sph_per(cen_snap_ray,0.05*r_par,res)
 chi_r_v=chi_r.vector().get_local()
 return([r_par,chi_r_v])

#if geo_p=='centre':
#ray_snap_cen=0.25
#csr_list=[[0.5,0.5,0.05*k] for k in range(1,1+Nsnap)]
#c_par : paramètre scalaire pour la position du centre
def snap_cen(c_par):
 cen_snap_ray=csr_list[c_par-1]
 chi_c=snapshot_sph_per(cen_snap_ray,ray_snap_cen,res)
 chi_c_v=chi_c.vector().get_local()
 return([c_par,chi_c_v])

# Calcul des snapshots, sous forme vectorielle


if parallelize:
 # Génération parallèle des snapshots
 pool=multiprocessing.Pool(processes=8)
 if geo_p=='rayon':
  list_chi_n_v=pool.map(snap_ray,(n for n in range(1,1+Nsnap)))
 elif geo_p=='centre':
  list_chi_n_v=pool.map(snap_cen,(n for n in range(1,1+Nsnap)))
 ## enregistrement des données dans une liste
else:
 # Génération séquentielle des snapshots : à compléter
 list_chi_n_v=[]
 for n in range(1,1+Nsnap):
  if geo_p=='rayon':
   chi_n_v=snap_ray(n*0.05)
  elif geo_p=='centre':
   chi_n_v=snap_cen(n*0.05)
  list_chi_n_v.append(chi_n_v)

# Construction de la liste des snapshots vectorisés : cas d'un paramètre géométrique définissant un ordre - lien avec la porosité ; ou non.

list_chi_v=[]
if geo_p=='rayon' or config=='compl':
 for n in range(1,1+Nsnap):
  for i in range(0,Nsnap):
   if list_chi_n_v[i][0]==n:
    chi_n_v=list_chi_n_v[i][1]
    list_chi_v.append(chi_n_v)
else:
 for i in range(0,Nsnap):
  chi_n_v=list_chi_n_v[i][1]
  list_chi_v.append(chi_n_v)

# Liste des snapshots : sauvegarde, on précise l'identité de la machine qui a effectué le calcul

l_name='Lchi_'+str(Nsnap)+'_'+config+'_'+geo_p+'_'+ordo+'_'+computer

# sauvegarde de la liste des solutions indexées calculées avec la méthode des éléments finis
with sh.open(repertoire_parent+l_name) as l_sto:
    l_sto["maliste"] = list_chi_v

# Matrice des snapshots : plus tard, voir l'étape II

# Exploitation des solution du problème aux éléments finis

for n in range(1,1+Nsnap):
 # Extraction du snapshot de rang n
 chi_n_v=list_chi_v[n-1]
 # On crée un maillage pour réécrire les snapshots sous la forme de fonctions
 if config=='sph_un':
  if geo_p=='ray':
   cen=cen_snap_ray
   r=n*0.05
   if typ_msh=='gms':
    print("maillages_per/3D/cubesphere_periodique_triangle_"+str(int(round(100*r,2)))+".xml")
    mesh=Mesh("maillages_per/3D/cubesphere_periodique_triangle_"+str(int(round(100*r,2)))+".xml")
   else:
    mesh=creer_maill_sph(cen,r,res)
  elif geo_p=='cen':
   r=ray_snap_cen
   mesh=creer_maill_sph(csr_list[n-1],r,res)
 elif config=='cylindre unique':
  if geo_p=='rayon':
   top=top_snap_ray
   r=n*0.05
   mesh=creer_maill_cyl(top,r,res)
  elif geo_p=='axe':
   r=ray_snap_ax
   mesh=creer_maill_cyl(acr_list[n-1],r,res)
 else:
  if geo_p=='rayon de la sphère variable':
   r=0
  elif geo_p=='rayon du cylindre variable':
   r=0
 V_n=VectorFunctionSpace(mesh, 'P', 2, constrained_domain=PeriodicBoundary())
 # On restitue la forme fonctionnelle du snapshot courant
 chi_n=Function(V_n)
 chi_n.vector().set_local(chi_n_v)
 # Représentation graphique
 plot(chi_n, linewidth=0.55)
 if fig_todo=='aff':
  plt.show()
 else:
  plt.savefig("Figures3D/sol_"+str(n)+"_sur"+str(Nsnap)+config+'_'+geo_p+".png")
 plt.close()
 # Affichage des valeurs et erreurs de la solution périodique, quelle que soit la configuration
 #err_per_ind_01(chi_n,cen,r,npas_err)
 err_per_gr(cen,r,chi_n,npas_err,fig_todo)
 # Tenseur de diffusion homogénéisé
 ## Intégrale de khi sur le domaine fluide
 T_chi=np.zeros((3,3))
 for k in range(0,3):
  for l in range(0,3):
   T_chi[k,l]=assemble(grad(chi_n)[k,l]*dx)
 ## Intégrale de l'identité sur le domaine fluide
 if config=='sphère unique':
  D=(1-4/3*pi*r**3)*np.eye(3)
 elif config=='cylindre unique':
  D=(1-pi*r**2)*np.eye(3)
 else :
  D=(1-4/3*pi*r_s**3-pi*r_c**2)*np.eye(3)
 ## Calcul et affichage du tenseur Dhom
 Dhom_k=D_k*(D+T_chi.T)
 #print(('Tenseur Dhom_k',Dhom_k))
 print('Coefficient Dhom_k11, snapshot '+str(n)+", "+config+', '+geo_p+" variable :",Dhom_k[0,0])
#
