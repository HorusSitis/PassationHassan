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

for cen in [[0.5,0.5],[0.0,0.0],[0.5,0.0],[0.0,0.5]]:
 mesh_c_r=creer_maill_circ(cen,r,res)
 #
 if fig_todo=='aff':
  plot(mesh_c_r)
  plt.show()
  plt.close()
 #
 elif fig_todo=='save':
  plot(mesh_c_r)
  plt.savefig("Figures2D/mesh_r_per"+str(int(round(100*cen[0],2)))+str(int(round(100*cen[1],2)))+str(int(round(100*r,2)))+".png")
  plt.close()

## Boucle pour la création des snapshots, avec un paramètre pouvant être le rayon d'une inclusion circulaire, ou l'emplacement de son centre ##
# Calcule aussi le tenseur de diffusion homogénéisé #

kfic=1
for i in range(7,1+Nsnap):#[0.111,0.211,0.316,0.423]:#,0.49]:#attention le rayon d'un cercle doit être non nul
 r=i*0.05
 c_x=0.5
 c_y=0.5
 chi_i=snapshot_circ_per([c_x,c_y],r,res)
 # Stockage des résultats avec un format hdf5
 ##LE.ecriture_champ_hdf5(kh_file,KH_SAVE,khi_i,kfic,file_rayon_ecriture,r,[c_x,c_y],res)
 print('Rayon :',r)
 print("Centre : "+str(c_x)+"_"+str(c_y))
 #fig_chi([c_x,c_y],r,chi_i,'aff')
 #fig_dchi([c_x,c_y],r,-grad(chi_i),fig_todo)
 err_per_gr([c_x,c_y],r,chi_i,50,fig_todo)
 #err_per_ind_01(chi_i,20)
 ##
 # Tenseur de diffusion homogénéisé
 ## Intégrale de khi sur le domaine fluide
 T_chi=array([[0.,0.],[0.,0.]])
 for k in range(0,2):
  for l in range(0,2):
   T_chi[k,l]=assemble(grad(chi_i)[k,l]*dx)
 ## Intégrale de l'identité sur le domaine fluide
 D=(1-pi*r**2)*np.eye(2)
 ## Calcul et affichage du tenseur Dhom
 Dhom_k=D_k*(D+T_chi.T)
 #print(('Tenseur D_hom',Dhom))
 print(Dhom_k[0,0])
 # Stockage
 ## ...
