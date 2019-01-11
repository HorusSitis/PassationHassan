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
mesh_fixe=generate_mesh(domaine_fixe,res_fixe)
V_fixe=VectorFunctionSpace(mesh_fixe, "P", 2, constrained_domain=PeriodicBoundary())

if fixe_aff:
 #représentation graphique du maillage
 plot(mesh_fixe)
 plt.show()
 plt.close()
else:
 #sauvegarde de la figure
 plot(mesh_fixe)
 plt.savefig("Figures3D/mesh_fixe.png")
 plt.close()

###-------------------- Commandes pour l'écriture de fichiers, à déplacer dans le script éventuellement --------------------###

#if [c_x,c_y]==[0.5,0.5]:
# suffixe="inc_centre/"
#elif [c_x,c_y]==[0.0,0.0]:
# suffixe="coins/"

#repertoire_final=repertoire_parent+suffixe
#File(repertoire_parent+"mesh_ulaire.xml.gz") << mesh_fixe
#ch_file,KH_SAVE=creation_fichier_pourecriture_champ_hdf5(repertoire_final,mesh_fixe)
#file_rayon_ecriture = open("%s/rayon_ecriture.txt" %(repertoire_final), "w")#kfic=1
# Famille de cellules élémentaires : 8 clichés, inclusion circulaire, paramétrée par le rayon du cercle
# Exemples de maillages raffiniés autour d'une inclusion périodique

#r=0.25
for cen in [[0.5,0.5,0.5]]:#,[0.0,0.0,0.0],[0.5,0.0,0.5],[0.0,0.5,0.0]]:
 mesh_s_r=creer_maill_sph(cen,0.25,res)
 #
 if fig_todo=='aff':
  plot(mesh_s_r)
  plt.show()
  plt.close()
 #
 elif fig_todo=='save':
  plot(mesh_s_r)
  plt.savefig("Figures3D/mesh_r_per"+str(int(round(100*cen[0],2)))+str(int(round(100*cen[1],2)))+str(int(round(100*cen[2],2)))+str(int(round(100*r,2)))+".png")
  plt.close()

## Boucle pour la création des snapshots, avec un paramètre pouvant être le rayon d'une inclusion circulaire, ou l'emplacement de son centre 

# Pour avoir des fonctions "top-level" à paralléliser

##res=12

#if geo_p=='rayon':
cen_snap_ray=[0.5,0.5,0.5]
def snap_ray(r_par):
 chi_r=snapshot_sph_per(cen_snap_ray,0.05*r_par,res)
 chi_r_v=chi_r.vector().get_local()
 return([r_par,chi_r_v])

#if geo_p=='centre':
ray_snap_cen=0.25
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
  list_chi_v=pool.map(snap_ray,(n for n in range(1,1+Nsnap)))
 elif geo_p=='centre':
  list_chi_v=pool.map(snap_cen,(n for n in range(1,1+Nsnap)))
 ## enregistrement des données dans une liste
else:
 # Génération séquentielle des snapshots : à compléter
 list_chi_v=[]
 for n in range(1,1+Nsnap):
  if geo_p=='rayon':
   chi_n_v=snap_ray(n*0.05)
  elif geo_p=='centre':
   chi_n_v=snap_cen(n*0.05)
  list_chi_v.append(chi_n_v)

# Liste des snapshots ...

# Stockage
## LE.ecriture_champ_hdf5(kh_file,KH_SAVE,khi_i,kfic,file_rayon_ecriture,r,[c_x,c_y],res)
## ...

#sys.exit()

# Exploitation des solution du problème aux éléments finis

for n in range(1,1+Nsnap):
 # Extraction du snapshot de rang n
 for i in range(0,Nsnap):
  if list_chi_v[i][0]==n:
   chi_n_v=list_chi_v[i][1]
 ## On crée un maillage pour réécrire les snapŝhots sous la forme de fonctions
 if geo_p=='rayon':
  cen=cen_snap_ray
  mesh=creer_maill_sph(cen,n*0.05,res)
 elif geo_p=='centre':
  r=ray_snap_cen
  mesh=creer_maill_sph(csr_list[n-1],r,res)
 V_n=VectorFunctionSpace(mesh, 'P', 2, constrained_domain=PeriodicBoundary())
 chi_n=Function(V_n)
 chi_n.vector().set_local(chi_n_v)
 # Affichage des valeurs et erreurs de la solution périodique
 ### procédures à écrire dans DDD_fun_obj.py
 ##
 # Tenseur de diffusion homogénéisé
 ## Intégrale de khi sur le domaine fluide
 T_chi=np.zeros((3,3))
 for k in range(0,3):
  for l in range(0,3):
   T_chi[k,l]=assemble(grad(chi_n)[k,l]*dx)
 ## Intégrale de l'identité sur le domaine fluide
 D=(1-4/3*pi*r**3)*np.eye(3)
 ## Calcul et affichage du tenseur Dhom
 Dhom_k=D_k*(D+T_chi.T)
 #print(('Tenseur Dhom_k',Dhom_k))
 print('Coefficient Dhom_k11, snapshot '+str(n)+", "+geo_p+" variable :",Dhom_k[0,0])
