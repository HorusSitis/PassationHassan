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

#if [c_x,c_y]==[0.5,0.5]:
# suffixe="inc_centre/"
#elif [c_x,c_y]==[0.0,0.0]:
# suffixe="coins/"

#repertoire_final=repertoire_parent+suffixe

#File(repertoire_parent+"mesh_circulaire.xml.gz") << mesh_fixe

#ch_file,KH_SAVE=creation_fichier_pourecriture_champ_hdf5(repertoire_final,mesh_fixe)

#file_rayon_ecriture = open("%s/rayon_ecriture.txt" %(repertoire_final), "w")
#kfic=1

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

#if geo_p=='rayon':
#cen_snap_ray=[0.5,0.5]

#if geo_p=='centre':
#ray_snap_cen=0.25
#csr_list=[[0.05*k,0.5]] for k in range(1,1+Nsnap)]
#c_par : paramètre scalaire pour la position du centre

for n in range(1,1+Nsnap):#[0.111,0.211,0.316,0.423]:#,0.49]:#attention le rayon d'un cercle doit être non nul
 if config=='cercle unique':
  if geo_p=='rayon':
   chi_n=snapshot_circ_per(cen_snap_ray,0.05*n,res) 
   c_x=cen_snap_ray[0]
   c_y=cen_snap_ray[1]
  elif geo_p=='centre':
   cen_snap_ray=csr_list[n-1]
   c_x=cen_snap_ray[0]
   c_y=cen_snap_ray[1]
   chi_n=snapshot_circ_per(cen_snap_ray,ray_snap_cen,res)
 #elif config=='cercles en alignés en diagonale':
 #elif config=='cercles alignés horizontalement':
 # figures et erreurs
 fig_chi([c_x,c_y],0.05*n,chi_n,fig_todo)
 #fig_dchi([c_x,c_y],r,-grad(chi_n),fig_todo)
 #err_per_gr([c_x,c_y],r,chi_n,50,fig_todo)
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
 #print(('Tenseur Dhom_k',Dhom_k))
 print('Coefficient Dhom_k11, snapshot '+str(n)+", "+config+', '+geo_p+" variable :",Dhom_k[0,0])
 # Stockage
 ## ...
