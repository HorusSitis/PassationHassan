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
 # maillage du domaine fixe avec gmsh
 if dom_fixe=='':
  mesh_fixe=Mesh("maillages_per/3D/cubesphere_periodique_triangle.xml")
 elif dom_fixe=="am":
  mesh_fixe=Mesh("maillages_per/3D/cubesphere_periodique_triangle"+"_am"+"_sur"+str(res)+"_fixe.xml")
 elif dom_fixe=='0001':
  mesh_fixe=Mesh("maillages_per/3D/cubesphere_periodique_triangle_sur"+str(res)+"_0001fixe.xml")
 elif dom_fixe=='0000':
  mesh_fixe=Mesh("maillages_per/3D/cubesphere_periodique_triangle_sur"+str(res)+"_0000fixe.xml")
else:
 mesh_fixe=generate_mesh(domaine_fixe,res_fixe)

V_fixe=VectorFunctionSpace(mesh_fixe, "P", 2, constrained_domain=PeriodicBoundary())

#sys.exit("chargement du maillage fixe terminé")
## Boucle pour la création des snapshots, avec un paramètre pouvant être le rayon d'une inclusion circulaire, ou l'emplacement de son centre 

# Pour avoir des fonctions "top-level" à paralléliser

## Sphère unique

#if geo_p=='ray':
cen_snap_ray=[0.5,0.5,0.5]
def snap_sph_ray(r_par):
 chi_r=snapshot_sph_per(cen_snap_ray,0.05*r_par,res_gmsh)
 chi_r_v=chi_r.vector().get_local()
 return([r_par,chi_r_v])

#if geo_p=='cen':
#ray_snap_cen=0.25
#csr_list=[[0.5,0.5,0.05*k] for k in range(1,1+Nsnap)]
#c_par : paramètre scalaire pour la position du centre
def snap_sph_cen(c_par):
 cen_snap_ray=csr_list[c_par-1]
 chi_c=snapshot_sph_per(cen_snap_ray,ray_snap_cen,res_gmsh)
 chi_c_v=chi_c.vector().get_local()
 return([c_par,chi_c_v])

## Cylindre unique

#if geo_p=='ray':
axe_snap_ray=[0.5,0.5]#,0.5]
def snap_cyl_ray(r_par):
 chi_r=snapshot_sph_per(axe_snap_ray,0.05*r_par)
 chi_r_v=chi_r.vector().get_local()
 return([r_par,chi_r_v])

#if geo_p=='axe':
#ray_snap_axe=0.25
asr_list=[[0.5,0.3+0.05*k] for k in range(1,1+Nsnap)]
#c_par : paramètre scalaire pour la position du centre
def snap_cyl_axe(c_par):
 cen_snap_ray=csr_list[c_par-1]
 chi_c=snapshot_sph_per(axe_snap_ray,ray_snap_axe,res)
 chi_c_v=chi_c.vector().get_local()
 return([c_par,chi_c_v])

# ------------------------- Snapshots, conditionnellement ------------------------- #

if not snap_done:
 # Calcul des snapshots, sous forme vectorielle
 ##if parallelize: Calcul parallèle oblligatoire
 # Génération parallèle des snapshots
 pool=multiprocessing.Pool(processes=8)
 if config=='sph_un':
  if geo_p=='ray':
   list_chi_n_v=pool.map(snap_sph_ray,(n for n in range(1,1+Nsnap)))
  elif geo_p=='cen':
   list_chi_n_v=pool.map(snap_sph_cen,(n for n in range(1,1+Nsnap)))
 elif config=='cyl_un':
  if geo_p=='ray':
   list_chi_n_v=pool.map(snap_cyl_ray,(n for n in range(1,1+Nsnap)))
  elif geo_p=='axe':
   list_chi_n_v=pool.map(snap_cyl_axe,(n for n in range(1,1+Nsnap)))
 ## enregistrement des données dans une liste
 # Construction de la liste des snapshots vectorisés : cas d'un paramètre géométrique définissant un ordre - lien avec la porosité ; ou non.
 list_chi_v=[]
 if geo_p=='ray' or config=='compl':
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
 l_name='Lchi_'+str(Nsnap)+'_'+config+'_'+geo_p+'_'+"sur"+str(res)+'_'+ordo+'_'+computer
 # sauvegarde de la liste des solutions indexées calculées avec la méthode des éléments finis
 with sh.open(repertoire_parent+l_name) as l_sto:
  l_sto["maliste"] = list_chi_v
 # Matrice des snapshots : plus tard, voir l'étape II
else :
 l_name='Lchi_'+str(Nsnap)+'_'+config+'_'+geo_p+'_'+"sur"+str(res)+'_'+ordo+'_'+computer
 with sh.open(repertoire_parent+l_name) as l_loa:
  list_chi_v = l_loa["maliste"]

# --------------------------------------------------------------------------------- #

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
    print("maillages_per/3D/cubesphere_periodique_triangle_"+str(int(round(100*r,2)))+"sur"+str(res)+".xml")
    mesh=Mesh("maillages_per/3D/cubesphere_periodique_triangle_"+str(int(round(100*r,2)))+"sur"+str(res)+".xml")
   else:
    mesh=creer_maill_sph(cen,r,res)
  elif geo_p=='cen':
   r=ray_snap_cen
   mesh=creer_maill_sph(csr_list[n-1],r,res)
 elif config=='cyl_un':
  if geo_p=='ray':
   top=top_snap_ray
   r=n*0.05
   mesh=creer_maill_cyl(top,r,res)
  elif geo_p=='axe':
   r=ray_snap_ax
   mesh=creer_maill_cyl(acr_list[n-1],r,res)
 else:
  if geo_p=='ray_sph':
   r=0
  elif geo_p=='ray_cyl':
   r=0
 V_n=VectorFunctionSpace(mesh, 'P', 2, constrained_domain=PeriodicBoundary())
 # On restitue la forme fonctionnelle du snapshot courant
 chi_n=Function(V_n)
 chi_n.vector().set_local(chi_n_v)
 # Représentation graphique
 plot(chi_n, linewidth=0.27)#35)
 if n==1:
  plt.title("R = 0,05", fontsize=40)
 else:
  plt.title("R = 0,"+str(int(round(100*r,2))),fontsize=40)
 if fig_todo=='aff':
  plt.show()
 else:
  plt.savefig("Figures3D/sol_"+str(n)+"_sur"+str(Nsnap)+config+'_'+geo_p+"res"+str(res)+".png")
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
 if config=='sph_un':
  D=(1-4/3*pi*r**3)*np.eye(3)
 elif config=='cyl_un':
  D=(1-pi*r**2)*np.eye(3)
 else :
  D=(1-4/3*pi*r_s**3-pi*r_c**2)*np.eye(3)
 ## Calcul et affichage du tenseur Dhom
 Dhom_k=D_k*(D+T_chi.T)
 #print(('Tenseur Dhom_k',Dhom_k))
 print("Noeuds",V_n.dim())
 print('Coefficient Dhom_k11EF, snapshot '+str(n)+", "+conf_mess+', '+geo_mess+" :",Dhom_k[0,0])
#
