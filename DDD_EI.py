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

import time

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
axe_snap_ray=[0.5,0.5]
def snap_cyl_ray(r_par):
 chi_r=snapshot_cyl_per(axe_snap_ray,0.05*r_par,res_gmsh)
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

def link(rho):
 por=1-(4/3*pi*0.35**3+pi*0.15**2)
 ray=1
 return(ray)

r_s_0=0
r_c_0=0
r_v_0=0

if config=='2sph':
 if deb==1:
  r_v_0=0.15
 else:
  r_v_0=0.2
elif config=='cylsph':
 if geo_p=='ray_cyl':
  if deb==1:
   r_s_0=0.15
  else:
   r_s_0=0.35
 elif geo_p=='ray_sph':
  if deb==1:
   r_c_0=0.15
  else:
   r_c_0=0.25
 elif geo_p=='ray_linked':
  print('aucun rayon fixé')

def snap_compl_ray(rho_par):
 ## deux sphères ##
 if geo_p=='ray':
  rho=0.05*rho_par
  chi_compl=snapshot_compl_per(rho,r_v_0,config,res_gmsh)
 ## un cylindre et une sphère ##
 elif geo_p=='ray_sph':
  rho=0.05#*rho_par
  chi_compl=snapshot_compl_per(rho,r_c_0,config,res_gmsh)
 elif geo_p=='ray_cyl':
  rho=0.05*rho_par
  chi_compl=snapshot_compl_per(r_s_0,rho,config,res_gmsh)
 elif geo_p=='ray_linked':
  rho=0.15+0.02*(rho_par-1)
  ray_link=link(rho)#((1/3)*(0.35**3-rho**3)+0.25**2)**(0.5)
  chi_compl=snapshot_compl_per(rho,ray_link,config,res_gmsh)
 ## on vectorise la fonction calculée par MEF ##
 chi_compl_v=chi_compl.vector().get_local()
 ## on renvoie un vecteur étiqueté, utilisable avec l'option 'par8' ##
 return([rho_par,chi_compl_v])






# ------------------------- Snapshots, conditionnellement ------------------------- #

if not snap_done:
 # -------- Calcul des snapshots, sous forme vectorielle, avec des étiquettes -------- #
 ### Génération parallèle des snapshots ###
 if gen_snap=='par8':
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
  else:
   list_chi_n_v=pool.map(snap_compl_ray,(n for n in range(1,1+Nsnap)))
 #elif config=='2sph':
 #elif config=='cylsph':
 # if geo_p=='ray_sph':
 # elif geo_p=='ray_cyl':
 # elif geo_p=='ray_linked':
 ### Génération séquentielle des snapshots, pour des tests de la méthode des éléments finis ###
 elif gen_snap=='seq':
  start=time.time()
  list_chi_n_v=[]
  for n in range(deb,deb+Nsnap):
   print(n)
   if config=='sph_un':
    if geo_p=='ray':
     list_chi_n_v.append(snap_sph_ray(n))
    elif geo_p=='cen':
     list_chi_n_v.append(snap_sph_cen(n))
   elif config=='cyl_un':
    if geo_p=='ray':
     list_chi_n_v.append(snap_cyl_ray(n))
    elif geo_p=='axe':
     list_chi_n_v.append(snap_cyl_axe(n))
   else:
    list_chi_n_v.append(snap_compl_ray(n))
  end=time.time()
  print('temps EF : ',end-start,' secondes')
 ### Génération parallèle pour chaque snapshot, pour de gros maillages ###
 elif gen_snap=='seq_par':
  list_chi_n_v=[]
 # -------- enregistrement des fonctions vectorisées dans une liste -------- #
 # Construction de la liste des snapshots vectorisés : cas d'un paramètre géométrique définissant un ordre - lien avec la porosité ; ou non.
 list_chi_v=[]
 if deb!=1:
  list_chi_v.append(list_chi_n_v[0][1])
 elif geo_p=='ray' or config=='compl':
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
res=res_gmsh
for n in range(deb,deb+Nsnap):
 # Extraction du snapshot de rang n
 chi_n_v=list_chi_v[n-deb]
 # On crée un maillage pour réécrire les snapshots sous la forme de fonctions
 if config=='sph_un':
  if geo_p=='ray':
   cen=cen_snap_ray
   r=n*0.05
   if typ_msh=='gms':
    mesh_name="cubesphere_periodique_triangle_"+str(int(round(100*r,2)))+"sur"+str(res)
   else:
    mesh=creer_maill_sph(cen,r,res)
  elif geo_p=='cen':
   r=ray_snap_cen
   mesh=creer_maill_sph(csr_list[n-1],r,res)
 elif config=='cyl_un':## avec gmsh
  if geo_p=='ray':
   r=n*0.05
   mesh_name="cubecylindre_periodique_triangle_"+str(int(round(100*r,2)))+"sur"+str(res)
  #elif geo_p=='axe':
   #r=ray_snap_ax
   #mesh=creer_maill_cyl(acr_list[n-1],r,res)
 elif config=='2sph':
  r=n*0.05
  r_s=r
  r_v=r_v_0
  mesh_name="cube"+config+"_periodique_triangle_"+str(int(round(100*r_s,2)))+str(int(round(100*r_v_0,2)))+"sur"+str(res)
 elif config=='cylsph':
  r=n*0.05
  if geo_p=='ray_cyl':
   r_c=r
   r_s=r_s_0
  elif geo_p=='ray_sph':
   r_c=r_c_0
   r_s=r
  mesh_name="cube"+config+"_periodique_triangle_"+str(int(round(100*r_c,2)))+str(int(round(100*r_s,2)))+"sur"+str(res)
 print("maillages_per/3D/"+mesh_name+".xml")
 mesh=Mesh("maillages_per/3D/"+mesh_name+".xml")
 V_n=VectorFunctionSpace(mesh, 'P', 2, constrained_domain=PeriodicBoundary())
 # On restitue la forme fonctionnelle du snapshot courant
 chi_n=Function(V_n)
 chi_n.vector().set_local(chi_n_v)
 # Représentation graphique
 plot(chi_n, linewidth=0.27)#35)
 plt.tight_layout(pad=0)
 if r<0.1:
  plt.title("Rho = 0,05", fontsize=40)
 elif deb==1:
  plt.title("Rho = 0,"+str(int(round(100*r,2))),fontsize=40)
 if fig_todo=='aff':
  plt.show()
 else:
  plt.savefig("Figures3D/sol_"+str(n)+"_sur"+str(Nsnap)+config+'_'+geo_p+"res"+str(res)+".png")
 plt.close()
 # Affichage des valeurs et erreurs de la solution périodique, quelle que soit la configuration
 #err_per_ind_01(chi_n,cen,r,npas_err)
 if config=='cyl_un' and geo_p=='ray':
  cen=[0.5,0.,0.5]# on triche un peu : on prend une face prévée d'une demie-sphère au lieu d'une face privée du disque frontal du cylindre
 if config=='2sph':
  err_per_gr_compl(config,r_v,chi_n,npas_err,fig_todo)
 elif config=='cylsph':
  err_per_gr_compl(config,r_c,chi_n,npas_err,fig_todo)
 else:
  err_per_gr(cen,r,chi_n,npas_err,fig_todo)
 # Tenseur de diffusion homogénéisé
 ## Intégrale de khi sur le domaine fluide
 T_chi=np.zeros((3,3))
 for k in range(0,3):
  for l in range(0,3):
   T_chi[k,l]=assemble(grad(chi_n)[k,l]*dx)
 #print(T_chi)
 ## Intégrale de l'identité sur le domaine fluide
 if config=='sph_un':
  por=1-4/3*pi*r**3
 elif config=='cyl_un':
  por=1-pi*r**2
 elif config=='2sph':
  por=1-4/3*pi*(r_s**3+r_v**3)
 elif config=='cylsph' :
  por=1-4/3*pi*r_s**3-pi*r_c**2
 D=por*np.eye(3)
 ## Calcul et affichage du tenseur Dhom
 Dhom_k=D_k*(D+T_chi.T)
 #print(('Tenseur Dhom_k',Dhom_k))
 print("Noeuds",V_n.dim())
 print("Porosité :",por)
 print('Coefficient Dhom_k11EF, snapshot '+str(n)+", "+conf_mess+', '+geo_mess+" :",Dhom_k[0,0])
 integ=assemble(chi_n[1]*dx)
 print('Valeur moyenne : ',integ)
#
