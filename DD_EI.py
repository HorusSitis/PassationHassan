#######################################################################################################################################
## Etape I : réalisation des clichés, avec la méthode des éléments finis. Calcul du tenseur d'homogénéisation. Stockage dans snap2D/ ##
#######################################################################################################################################

### ------------ Reproduire éventuellement pour des étapes ultérieures. Laisser seulement dans DD_fun_obj ? ------------ ###

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

#res_fixe=30#résolution du maillage sans obstacle

domaine_fixe=Rectangle(Point(xinf,yinf),Point(xsup,ysup))
mesh_fixe=generate_mesh(domaine_fixe,res_fixe)
V_fixe=VectorFunctionSpace(mesh_fixe, "P", 3, form_degree=1, constrained_domain=PeriodicBoundary())

if fixe_aff:
 #représentation graphique du maillage
 plot(mesh_fixe)
 plt.show()
 plt.close()
else:
 #sauvegarde de la figure
 plot(mesh_fixe)
 plt.tight_layout()
 plt.savefig("Figures2D/mesh_fixe.png")
 plt.close()

###-------------------- Commandes pour l'écriture de fichiers, à déplacer dans le script éventuellement --------------------###



#repertoire_final=repertoire_parent+suffixe
#File(repertoire_parent+"mesh_circulaire.xml.gz") << mesh_fixe
#ch_file,KH_SAVE=creation_fichier_pourecriture_champ_hdf5(repertoire_final,mesh_fixe)
#file_rayon_ecriture = open("%s/rayon_ecriture.txt" %(repertoire_final), "w")#kfic=1
# Famille de cellules élémentaires : 8 clichés, inclusion circulaire, paramétrée par le rayon du cercle
# Exemples de maillages raffiniés autour d'une inclusion périodique

r=0.25

for cen in []:#[[0.5,0.5],[0.0,0.0],[0.5,0.0],[0.0,0.5]]:
 mesh_c_r=creer_maill_circ(cen,r,res)
 #
 if fig_todo=='aff':
  plot(mesh_c_r)
  plt.show()
  plt.close()
 #
 elif fig_todo=='save':
  plot(mesh_c_r)
  plt.tight_layout()
  plt.savefig("Figures2D/mesh_r_per"+str(int(round(100*cen[0],2)))+str(int(round(100*cen[1],2)))+str(int(round(100*r,2)))+".png")
  plt.close()

#sys.exit()#---------------------------------------------------

## Boucle pour la création des snapshots, avec un paramètre pouvant être le rayon d'une inclusion circulaire, ou l'emplacement de son centre ##
# Calcule aussi le tenseur de diffusion homogénéisé #


## Cercle unique

#if geo_p=='ray':
cen_snap_ray=[0.5,0.5]#,0.5]
def snap_circ_ray(r_par):
 chi_r=snapshot_circ_per(cen_snap_ray,0.05*r_par,res)
 chi_r_v=chi_r.vector().get_local()
 return([r_par,chi_r_v])

#if geo_p=='cen':
#ray_snap_cen=0.25
#csr_list=[[0.5,0.3+0.05*k] for k in range(1,1+Nsnap)]
#c_par : paramètre scalaire pour la position du centre
def snap_circ_cen(c_par):
 cen_snap_ray=csr_list[c_par-1]
 chi_c=snapshot_circ_per(cen_snap_ray,ray_snap_cen,res)
 chi_c_v=chi_c.vector().get_local()
 return([c_par,chi_c_v])

# ------------------------- Snapshots, conditionnellement ------------------------- #

if not snap_done:
 # Calcul des snapshots, sous forme vectorielle
 ##if parallelize: Calcul parallèle oblligatoire
 # Génération parallèle des snapshots
 pool=multiprocessing.Pool(processes=8)
 if config=='cer_un':
  if geo_p=='ray':
   list_chi_n_v=pool.map(snap_circ_ray,(n for n in range(1,1+Nsnap)))
  elif geo_p=='cen':
   list_chi_n_v=pool.map(snap_circ_cen,(n for n in range(1,1+Nsnap)))
 #elif config=='cyl_un':
 # if geo_p=='ray':
 #  list_chi_n_v=pool.map(snap_cyl_ray,(n for n in range(1,1+Nsnap)))
 # elif geo_p=='axe':
 #  list_chi_n_v=pool.map(snap_cyl_axe,(n for n in range(1,1+Nsnap)))
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
 l_name='Lchi_'+str(Nsnap)+'_'+config+'_'+geo_p+'_'+ordo+'_'+computer
 # sauvegarde de la liste des solutions indexées calculées avec la méthode des éléments finis
 with sh.open(repertoire_parent+l_name) as l_sto:
  l_sto["maliste"] = list_chi_v
 # Matrice des snapshots : plus tard, voir l'étape II
else :
 l_name='Lchi_'+str(Nsnap)+'_'+config+'_'+geo_p+'_'+ordo+'_'+computer
 with sh.open(repertoire_parent+l_name) as l_loa:
  list_chi_v = l_loa["maliste"]

# --------------------------------------------------------------------------------- #

for n in range(1,1+Nsnap):
 r=0.05*n
 mesh=Mesh("maillages_per/2D/maillage_trou2d_"+str(int(round(100*r,2)))+".xml")
 V_n=VectorFunctionSpace(mesh, 'P', 2, constrained_domain=PeriodicBoundary())
 print(str(n),V_n.dim())





for n in range(1,1+Nsnap):#[0.111,0.211,0.316,0.423]:#,0.49]:#attention le rayon d'un cercle doit être non nul
 # Extraction du snapshot de rang n
 chi_n_v=list_chi_v[n-1]
 # On crée un maillage pour réécrire les snapshots sous la forme de fonctions
 if config=='cer_un':
  if geo_p=='ray':
   cen=cen_snap_ray
   r=n*0.05
   if typ_msh=='gms':
    #print("maillages_per/3D/cubesphere_periodique_triangle_"+str(int(round(100*r,2)))+".xml")
    mesh=Mesh("maillages_per/2D/maillage_trou2d_"+str(int(round(100*r,2)))+".xml")
   else:
    mesh=creer_maill_circ(cen,r,res)
  #elif geo_p=='cen':
 V_n=VectorFunctionSpace(mesh, 'P', 2, constrained_domain=PeriodicBoundary())
 # On restitue la forme fonctionnelle du snapshot courant
 chi_n=Function(V_n)
 print(V_n.dim())
 print(len(chi_n_v))
 sys.exit()#-----------------------------------------------------------
 chi_n.vector().set_local(chi_n_v)
 # Figures et erreurs
 plot(chi_n)
 #plot(grad(chi_n)[:,0]
 #plot(grad(chi_n)[:,1]
 if fig_todo=='aff':
  plt.show()
 #elif fig_todo=='save':
 plt.close()
 #fig_chi([c_x,c_y],0.05*n,chi_n,fig_todo)
 err_per_gr([c_x,c_y],r,chi_n,npas_err,fig_todo)
 #err_per_ind_01(chi_n,20)
 #sys.exit()#---------------------------------------------------
 ##
 # Tenseur de diffusion homogénéisé
 ## Intégrale de chi sur le domaine fluide
 T_chi=array([[0.,0.],[0.,0.]])
 for k in range(0,2):
  for l in range(0,2):
   T_chi[k,l]=assemble(grad(chi_n)[k,l]*dx)
 ## Intégrale de l'identité sur le domaine fluide
 D=(1-pi*r**2)*np.eye(2)
 ## Calcul et affichage du tenseur Dhom
 Dhom_k=D_k*(D+T_chi.T)
 print(('Tenseur Dhom_k',Dhom_k))
 print('Coefficient Dhom_k11, snapshot '+str(n)+", "+config_mess+', '+geo_mess+" variable :",Dhom_k[0,0])
