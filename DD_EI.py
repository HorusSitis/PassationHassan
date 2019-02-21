#######################################################################################################################################
## Etape I : réalisation des clichés, avec la méthode des éléments finis. Calcul du tenseur d'homogénéisation. Stockage dans snap2D/ ##
#######################################################################################################################################

### ------------ Reproduire éventuellement pour des étapes ultérieures. Laisser seulement dans DD_fun_obj ? ------------ ###

tol=1e-10

xinf=0.0
yinf=0.0
xsup=1.0
ysup=1.0

import time

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

## Domaine fixe sans gmsh
#res_fixe=30#résolution du maillage sans obstacle
#domaine_fixe=Rectangle(Point(xinf,yinf),Point(xsup,ysup))
#mesh_fixe=generate_mesh(domaine_fixe,res_fixe)

if typ_msh=='gms':
 res_fixe=res_gmsh
 if dom_fixe=="am":
  mesh_f_name="maillages_per/2D/maillage_fixe2D_am.xml"
 elif config=='compl':
  mesh_f_name="maillages_per/2D/maillage_trous2D_"+geo_p+"_fixe.xml"

mesh_fixe=Mesh(mesh_f_name)

V_fixe=VectorFunctionSpace(mesh_fixe, "P", VFS_degree, form_degree=1, constrained_domain=PeriodicBoundary())

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



## Boucle pour la création des snapshots, avec un paramètre pouvant être le rayon d'une inclusion circulaire, ou l'emplacement de son centre ##
# Calcule aussi le tenseur de diffusion homogénéisé #


## Cercle unique

#if geo_p=='ray':
#cen_snap_ray=[0.5,0.5]
def snap_circ_ray(r_par):
 if test_snap=='i_per':
  chi_r=snapshot_circ_per(cen_snap_ray,0.05*r_par,res)
 else:
  chi_r=snapshot_compl_per(geo_p,0.05*r_par,cen_snap_ray,mention,test_snap)
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

def snap_compl_ray(r_par):
 if geo_p=='diag':
  rho=0.05*r_par
 elif geo_p=='hor':
  rho=0.01+0.04*r_par
 chi_compl=snapshot_compl_per(geo_p,rho,cen_snap_ray,mention,test_snap)
 chi_compl_v=chi_compl.vector().get_local()
 return([r_par,chi_compl_v])

# ------------------------- Snapshots, conditionnellement ------------------------- #
#sys.exit("test pour l'homogénéisation périodique effectué")#---------------------------------------------------
if not snap_done:
 # Calcul des snapshots, sous forme vectorielle
 if gen_snap=='par8':
  # Génération parallèle des snapshots
  pool=multiprocessing.Pool(processes=8)
  if config=='cer_un':
   if geo_p=='ray':
    list_chi_n_v=pool.map(snap_circ_ray,(n for n in range(1,1+Nsnap)))
   elif geo_p=='cen':
    list_chi_n_v=pool.map(snap_circ_cen,(n for n in range(1,1+Nsnap)))
  elif config=='cer_un_som':
   list_chi_n_v=pool.map(snap_circ_ray,(n for n in range(1,1+Nsnap)))
  elif config=='compl':
   list_chi_n_v=pool.map(snap_compl_ray,(n for n in range(1,1+Nsnap)))
 elif gen_snap=='seq':
  start=time.time()
  list_chi_n_v=[]
  for n in range(deb,deb+Nsnap):
   if geo_p=='ray':
    list_chi_n_v.append(snap_circ_ray(n))
   elif geo_p=='cen':
    list_chi_n_v.append(snap_circ_cen(n))
   elif config=='cer_un_som':
    list_chi_n_v.append(snap_circ_ray(n))
   elif config=='compl':
    list_chi_n_v.append(snap_compl_ray(n))
 #elif parallelize=='seq_par':
  end=time.time()
  print('résolution EF : ',end-start,' secondes')
  sys.exit('test individuel de temps d éxécution terminé')
 ### utilisé en 3D, pour le calcul parallèle d'un snapshot individuel
 #
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
 l_name=test_snap+'Lchi_'+str(Nsnap)+'_'+config+'_'+geo_p+'_deg'+str(VFS_degree)+'_'+ordo+'_'+computer
 # sauvegarde de la liste des solutions indexées calculées avec la méthode des éléments finis
 with sh.open(repertoire_parent+l_name) as l_sto:
  l_sto["maliste"] = list_chi_v
 # Matrice des snapshots : plus tard, voir l'étape II
else :
 l_name=test_snap+'Lchi_'+str(Nsnap)+'_'+config+'_'+geo_p+'_deg'+str(VFS_degree)+'_'+ordo+'_'+computer
 with sh.open(repertoire_parent+l_name) as l_loa:
  list_chi_v = l_loa["maliste"]

#sys.exit("liste des solutions EF créée")
# --------------------------------------------------------------------------------- #

#mesh_fixe=Mesh("maillages_per/2D/maillage_fixe2d_am.xml")
#V_fixe=VectorFunctionSpace(mesh_fixe, 'P', 3, constrained_domain=PeriodicBoundary())

for n in range(1,1+Nsnap):#attention le rayon d'un cercle doit être non nul
 # Extraction du snapshot de rang n
 chi_n_v=list_chi_v[n-1]
 print(config)
 # On crée un maillage pour réécrire les snapshots sous la forme de fonctions
 mesh_directory="maillages_per/2D/"
 if config=='cer_un':
  if geo_p=='ray':
   cen=cen_snap_ray
   r=n*0.05
   if typ_msh=='gms':
    mesh_name="maillage_trou2D_"+str(int(round(100*r,2)))
    print(mesh_name)
    mesh=Mesh(mesh_directory+mesh_name+".xml")
    plot(mesh)
    #plt.title("Periodical mesh", fontsize=30)
    plt.tight_layout(pad=0)
    if fig_todo=='aff':
     plt.show()
    elif fig_todo=='save' and r==0.25:
     plt.savefig('Figures2D/maillage_gmsh_per_'+config+geo_p+'_ray'+str(int(round(100*r,2)))+'png')
    else:
     print('pfffrrrhhhh !!!')
    plt.close()
   else:
    mesh=creer_maill_circ(cen,r,res)
 elif config=='cer_un_som':
  r=n*0.05
  mention="_som"
  if typ_msh=='gms':
   mesh_name="maillage_trou2D"+mention+"_"+str(int(round(100*r,2)))
   print(mesh_name)
   mesh=Mesh(mesh_directory+mesh_name+".xml")
   plot(mesh)
   #plt.title("Periodical mesh", fontsize=30)
   plt.tight_layout(pad=0)
   if fig_todo=='aff':
    plt.show()
   elif fig_todo=='save' and r==0.25:
    plt.savefig('Figures2D/maillage_gmsh_per_'+config+geo_p+'_ray'+str(int(round(100*r,2)))+".png")
   else:
    print('pfffrrrhhhh !!!')
   plt.close()
  else:
   mesh=creer_maill_circ([c_x,c_y],r,res)
  #elif geo_p=='cen':
 else:
  if geo_p=='diag':
   rho=0.05*n
  elif geo_p=='hor':
   rho=0.01+0.04*n
  r=rho
  r_fixe=0.15
  mesh_name="maillage_trous2D_"+geo_p+"_"+str(int(round(100*rho,2)))
  print(mesh_name)
  mesh=Mesh(mesh_directory+mesh_name+".xml")
  plot(mesh)
  #plt.title("Periodical mesh", fontsize=30)
  plt.tight_layout(pad=0)
  if fig_todo=='aff':
   plt.show()
  elif fig_todo=='save' and r==0.25:
   plt.savefig('Figures2D/maillage_gmsh_per_'+config+geo_p+'_ray'+str(int(round(100*r,2)))+".png")
  else:
   print('pfffrrrhhhh !!!')
  plt.close()
 V_n=VectorFunctionSpace(mesh, 'P', VFS_degree, constrained_domain=PeriodicBoundary())
 # On restitue la forme fonctionnelle du snapshot courant
 chi_n=Function(V_n)
 print(V_n.dim())
 print(len(chi_n_v))
 #sys.exit("code déboggué")#---------------------------------------------------
 chi_n.vector().set_local(chi_n_v)
 # Figures et erreurs
 plot(chi_n)
 if n==1:
  plt.title("Rho = 0,05", fontsize=40)
 else:
  plt.title("Rho = 0,"+str(int(round(100*r,2))),fontsize=40)
 #plot(grad(chi_n)[:,0]
 #plot(grad(chi_n)[:,1]
 if fig_todo=='aff':
  plt.show()
 elif fig_todo=='save':
  plt.savefig("Figures2D/sol_"+str(n)+"_sur"+str(Nsnap)+config+'_'+geo_p+".png")
 plt.close()
 if config!='compl':
  #fig_chi([c_x,c_y],0.05*n,chi_n,fig_todo)
  err_per_gr(cen_snap_ray,r,chi_n,npas_err,fig_todo)
 else:
  #fig_chi([c_x,c_y],r_fixe,chi_n,fig_todo)
  err_per_gr(cen_snap_ray,r_fixe,chi_n,npas_err,fig_todo)
 #err_per_ind_01(chi_n,20)
 ##
 # Tenseur de diffusion homogénéisé
 ## Intégrale de chi sur le domaine fluide
 T_chi=np.zeros((2,2))
 for k in range(0,2):
  for l in range(0,2):
   T_chi[k,l]=assemble(grad(chi_n)[k,l]*dx)
 ## Intégrale de l'identité sur le domaine fluide
 if config!='compl':
  por=(1-pi*r**2)
 else:
  por=1-pi*(r**2+0.15**2)
 D=por*np.eye(2)
 ## Calcul et affichage du tenseur Dhom
 Dhom_k=D_k*(D+T_chi.T)
 #print(('Tenseur Dhom_k',Dhom_k))
 print("Porosité :", por)
 print('Coefficient Dhom_k11, snapshot '+str(n)+", "+conf_mess+', '+geo_mess+" :",Dhom_k[0,0])
