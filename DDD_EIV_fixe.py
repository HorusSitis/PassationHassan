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

if res_gmsh==10:
 lw=0.27
elif res_gmsh==20:
 lw=0.15
elif res_gmsh==50:
 lw=0.01

r_s_0=0.15
r_v_0=0.15
r_c_0=0.15

r_min=0.05

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

mesh_dir="maillages_per/3D/"

## inclusions simples ou rayons liés
if dom_fixe=="am":
 mesh_f_name=mesh_dir+"cube_periodique_triangle"+"_"+dom_fixe+"_sur"+str(res_gmsh)+"_fixe.xml"
## inclusions multiples, unique rayon variable
elif dom_fixe=="solid":
 mesh_fixe_prefix=mesh_dir+"cube"+config+"_periodique_triangle_"
 if config=='2sph':
  mesh_f_name=mesh_fixe_prefix+"fixe"+str(int(round(100*r_v_0,2)))+"sur"+str(res_gmsh)+".xml"
 elif config=='cylsph':
  ## rayon du cylindre aux arètes ou de la sphère centrale fixés à 0.15 ##
  if geo_p=='ray_sph':
   mesh_f_name=mesh_fixe_prefix+str(int(round(100*r_c_0,2)))+"fixe"+"sur"+str(res_gmsh)+".xml"
  elif geo_p=='ray_cyl':
   if fixe_comp=='cyl_sph':
    mesh_f_name=mesh_fixe_prefix+"fixe"+str(int(round(100*r_s_0,2)))+"sur"+str(res_gmsh)+".xml"
   elif fixe_comp=='sph_un':
    mesh_f_name=mesh_dir+"cubesphere_periodique_triangle_"+str(int(round(100*r_s_0,2)))+"sur"+str(res_gmsh)+".xml"
elif dom_fixe=="ray_min":
 fixe_comp=True#utilisation du domaine fixe avec annulation du rayon du cylindre dans le fichier général
 if config=='cylsph':
  if geo_p=='ray_sph':
   mesh_f_name=mesh_dir+"cube"+config+"_periodique_triangle_"+str(int(round(100*r_c_0,2)))+str(int(round(100*r_min,2)))+"sur"+str(res_gmsh)+".xml"
  elif geo_p=='ray_cyl':
   mesh_f_name=mesh_dir+"cube"+config+"_periodique_triangle_"+str(int(round(100*r_min,2)))+str(int(round(100*r_s_0,2)))+"sur"+str(res_gmsh)+".xml"

mesh_fixe=Mesh(mesh_f_name)

# fonctions test du domaine fixe

V_fixe=VectorFunctionSpace(mesh_fixe,'P',2,constrained_domain=PeriodicBoundary())

# Performances

import time

### ------------ Etapes reproduites : dépendances directes de Main3D ------------ ###

#r_nouv=0.22
nb_modes=N_mor

#mesh_dir="maillages_per/3D/"
if config=='sph_un':
 mesh_n_name=mesh_dir+"cubesphere_periodique_triangle_"+str(int(round(100*r_nouv,2)))+"sur"+str(res_gmsh)
elif config=='cyl_un':
 mesh_n_name=mesh_dir+"cubecylindre_periodique_triangle_"+str(int(round(100*r_nouv,2)))+"sur"+str(res_gmsh)
if config=='2sph':
 mesh_n_name=mesh_dir+"cube"+config+"_periodique_triangle_"+str(int(round(100*r_nouv,2)))+str(int(round(100*r_v_0,2)))+"sur"+str(res_gmsh)
elif config=='cylsph':
 if geo_p=='ray_sph':
  mesh_n_name=mesh_dir+"cube"+config+"_periodique_triangle_"+str(int(round(100*r_c_0,2)))+str(int(round(100*r_nouv,2)))+"sur"+str(res_gmsh)
 elif geo_p=='ray_cyl':
  mesh_n_name=mesh_dir+"cube"+config+"_periodique_triangle_"+str(int(round(100*r_nouv,2)))+str(int(round(100*r_s_0,2)))+"sur"+str(res_gmsh)
 #elif geo_p=='ray_linked':

### ------------ Etapes spécifiques à la construction du modèle réduit ------------ ###

nb_modes=N_mor

r=r_nouv
# Création du domaine d'intégration sur le maillage fixe : pour le calcul de Dhom
class DomPhysFluide(SubDomain):
 def inside(self, x, on_boundary):
  return True if ((x[0]-0.5)**2+(x[1]-0.5)**2+(x[2]-0.5)**2>=r**2) else False
# marquage du domaine
dom_courant=DomPhysFluide()
subdomains=MeshFunction('size_t',mesh_fixe,mesh_fixe.topology().dim())
subdomains.set_all(1)
dom_courant.mark(subdomains,12829)
dxf=Measure("dx", domain=mesh_fixe, subdomain_data=subdomains)

## Taille du maillage du domaine fixe ##

nb_noeuds_fixe=V_fixe.dim()

## Chargement de la base POD complète

phi_name='Phi'+dom_fixe+'_dim'+str(Nsnap)+'_'+config+'_'+geo_p+'_'+"res"+str(res)+'_'+ordo+'_'+computer

print(phi_name)

with sh.open(repertoire_parent+phi_name) as phi_loa:
 Phi_prime_v = phi_loa["maliste"]

## Création de la base POD tronquée, sous forme vectorielle

Phi_mor=Phi_prime_v[:,range(0,nb_modes)]

### ------------ Fin ------------ ###

#print("Maillage fixe : ",mesh_f_name)
#print(mesh_n_name)

mesh_nouv=Mesh(mesh_n_name+".xml")
V_nouv=VectorFunctionSpace(mesh_nouv, "P", 2, constrained_domain=PeriodicBoundary())




# --------------------- SE1 : projection de la base POD sur le nouveau domaine --------------------- #

## On initialise le temps de calcul ##

start_se1=time.time()

# Initialisation des fonctions POD

phi_fixe=Function(V_fixe)

# Raffinement du maillage : pour le moment, cas d'une inclusion unique centrée
def width(i):
 return crow/i**3

start_refi=time.time()
mesh_r_fixe=mesh_fixe
r=r_nouv

for i in range(1,1+Nrefine):
 print('raffinement',i)
 markers = MeshFunction("bool", mesh_fixe, mesh_fixe.topology().dim())
 markers.set_all(False)
 for c in cells(mesh_fixe):
  #for f in facets(c):
  # if ((f.midpoint()[0]-cen_snap_ray[0])**2+(f.midpoint()[1]-cen_snap_ray[1])**2+(f.midpoint()[2]-cen_snap_ray[2])**2<=(r*(1+width(i)))**2) and ((f.midpoint()[0]-cen_snap_ray[0])**2+(f.midpoint()[1]-cen_snap_ray[1])**2+(f.midpoint()[2]-cen_snap_ray[2])**2>=(r*(1-width(i)))**2):
  #  markers[c]=True
  #for v in vertices(c):
  # if ((v.point().x()-cen_snap_ray[0])**2+(v.point().y()-cen_snap_ray[1])**2+(v.point().z()-cen_snap_ray[2])**2>=(r*(1-width(i)))**2) and ((v.point().x()-cen_snap_ray[0])**2+(v.point().y()-cen_snap_ray[1])**2+(v.point().z()-cen_snap_ray[2])**2<=(r*(1+width(i)))**2):
  #  markers[c]=True
  list_sgn=[]
  for v in vertices(c):
   if (v.point().x()-cen_snap_ray[0])**2+(v.point().y()-cen_snap_ray[1])**2+(v.point().z()-cen_snap_ray[2])**2>r**2:
    list_sgn.append(1)
   elif (v.point().x()-cen_snap_ray[0])**2+(v.point().y()-cen_snap_ray[1])**2+(v.point().z()-cen_snap_ray[2])**2==r**2:
    list_sgn.append(0)
   else:
    list_sgn.append(-1)
  # on marque les cellules qui coupent la frontière du domaine fluide virtuel
  if list_sgn==[1,1,1,1] or list_sgn==[-1,-1,-1,-1]:
   markers[c]=False
  else:
   markers[c]=True
 mesh_r_fixe=refine(mesh_r_fixe, markers, redistribute=True)
print(list_sgn)
end=time.time()
tps_refi=end-start_refi

print('raffinement fait',tps_refi,'secondes')

if fig_todo=='aff':
 plot(mesh_r_fixe)
 plt.title("Maillage raffiné autour de la frontière physique, rayon "+str(int(round(100*r_nouv,2))))
 plt.show()
 plt.close()

V_r_fixe=VectorFunctionSpace(mesh_r_fixe, "P", 2, constrained_domain=PeriodicBoundary())

nb_noeuds_r_fixe=V_r_fixe.dim()
Phi_r_fixe_v=np.zeros((nb_noeuds_r_fixe,nb_modes))

print('Noeuds du domaine fixe :',nb_noeuds_fixe)
print('rho :',r_nouv)
print('Noeuds du maillage raffiné :',nb_noeuds_r_fixe)

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
   plt.savefig("Figures2D/phi_nouv_"+str(n+1)+"_"+config+'_'+geo_p+".png")
  plt.close()
  # on range le vecteur de POD interpolée dans la matrice Phi_nouv_v
  Phi_r_fixe_v[:,n]=phi_r_fixe.vector().get_local()
 # on mesure le temps d'interpolation
 end=time.time()
 tps_interp=end-start_interp

print('interpolation faite',tps_interp,'secondes')
## On enregistre et imprime le temps d'éxécution de SE1

end=time.time()

print('se1 faite',end-start_se1,'secondes')
sys.exit('plus de bug')#-------------------------------------


# --------------------- SE2 : résolution du modèle réduit --------------------- #

## On réinitialise le temps de calcul ##

start=time.time()




if config=='sph_un' or config=='cyl_un':
 Coeff=calc_Ab_simpl_3D_ninterpol(mesh_f_name,config,geo_p,r_nouv,Phi_prime_v,nb_modes)
else:
 if config=='2sph':
  r_cen=r_nouv
  r_per=r_v_0
 elif config=='cylsph':
  if geo_p=='ray_sph':
   r_cen=r_nouv
   r_per=r_c_0
  elif geo_p=='ray_cyl':
   r_cen=r_s_0
   r_per=r_nouv
 Coeff=calc_Ab_compl_3D_ninterpol(mesh_f_name,config,geo_p,r_cen,r_per,Phi_prime_v,nb_modes)
#sys.exit("solveur ROM éxécuté pour des inclusions multiples")
A=Coeff[0]
b=Coeff[1]

print('A :',A,'b :',b)
## On résoud le modèle réduit

a_nouv=np.linalg.solve(A.T,-b)

## On enregistre et imprime le temps d'éxécution de SE2

end=time.time()

print('se2 faite ',end-start,' secondes')
# --------------------- SE3 : calcul du nouveau champ de vecteurs, affichage --------------------- #

## On réinitialise le temps de calcul ##

start=time.time()

## On initialise et affiche le champ chi_nouv

chi_prime_nouv_v=np.dot(Phi_prime_v,a_nouv)
chi_prime_nouv=Function(V_fixe)
chi_prime_nouv.vector().set_local(chi_prime_nouv_v)

plot(chi_prime_nouv, linewidth=lw)
plt.title("Rho = 0,"+str(int(round(100*r_nouv,2))),fontsize=40)
if fig_todo=='aff':
 plt.show()
elif fig_todo=='save':
 plt.savefig("Figures3D/sol_romfixe"+str(int(round(100*r_nouv,2)))+"_sur"+str(Nsnap)+config+'_'+geo_p+"res"+str(res)+".png")
plt.close()

## Exploitation du champ ainsi obtenu
r=r_nouv
rho=r_nouv


### Affichage des valeurs et erreurs de la solution périodique, quelle que soit la configuration ###

if config=='sph_un' or config=='cyl_un':
 #err_per_ind_01(chi_n,cen,r,npas_err)
 err_per_gr(cen_snap_ray,0.05,chi_prime_nouv,npas_err,fig_todo)
elif config=='2sph':
 err_per_gr_compl(config,r_v_0,chi_prime_nouv,npas_err,fig_todo)
elif config=='cylsph':
 if geo_p=='ray_sph':
  err_per_gr_compl(config,r_c_0,chi_prime_nouv,npas_err,fig_todo)
 elif geo_p=='ray_cyl':
  err_per_gr_compl(config,r_nouv,chi_prime_nouv,npas_err,fig_todo)
 #elif geo_p=='ray_linked':


### Tenseur de diffusion homogénéisé ###

## Domaine fluide
if config=='sph_un' and geo_p=='ray':
 class DomPhysFluide(SubDomain):
  def inside(self, x, on_boundary):
   return True if (x[0]**2+x[1]**2+x[2]**2>=r**2) else False
elif config=='cyl_un' and geo_p=='ray':
 class DomPhysFluide(SubDomain):
  def inside(self, x, on_boundary):
   return True if (x[0]**2+x[2]**2>=r**2) else False
elif config=='2sph':
 class DomPhysFluide(SubDomain):
  def inside(self, x, on_boundary):
   return True if (x[0]**2+x[1]**2+x[2]**2>=r_cen**2) else False
elif config=='cylsph':
 class DomPhysFluide(SubDomain):
  def inside(self, x, on_boundary):
   return True if ((geo_p=='ray_sph' and (x[0]**2+x[1]**2+x[2]**2>=r_cen**2)) or (geo_p=='ray_cyl' and x[0]**2+x[2]**2>=r_per**2 and (1-x[0])**2+x[2]**2>=r_per**2 and x[0]**2+(1-x[2])**2>=r_per**2 and (1-x[0])**2+(1-x[2])**2>=r_per**2)) else False
dom_courant=DomPhysFluide()
subdomains=MeshFunction('size_t',mesh_fixe,mesh_fixe.topology().dim())
subdomains.set_all(1)
dom_courant.mark(subdomains,12829)
dxf=Measure("dx", domain=mesh_fixe, subdomain_data=subdomains)
## Intégrale de chi sur le domaine fluide
T_chi=np.zeros((3,3))
for k in range(0,3):
 for l in range(0,3):
  T_chi[k,l]=assemble(grad(chi_prime_nouv)[k,l]*dxf(12829))
## Intégrale de l'identité sur le domaine fluide
### Calcul de la porosité
if config=='sph_un':
 por=1-4/3*pi*r_nouv**3
elif config=='cyl_un':
 por=1-pi*r_nouv**2
elif config=='2sph':
 por=1-4/3*pi*(r_nouv**3+r_v_0**3)
elif config=='cylsph':
 if geo_p=='ray_sph':
  r_s=r_nouv
  r_c=r_c_0
 elif geo_p=='ray_cyl':
  r_s=r_s_0
  r_c=r_nouv
 #elif geo_p=='ray_linked':
 por=1-4/3*pi*r_s**3-pi*r_c**2
### Intégration du terme constant du coefficient d diffusion, sur le domaine fluide
D=por*np.eye(3)
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
if config=='sph_un':
 chi_nouv=snapshot_sph_per(cen_snap_ray,r_nouv,res,typ_sol)
elif config=='cyl_un':
 chi_nouv=snapshot_cyl_per(top_snap_ray,r_nouv,res,typ_sol)
elif config=='2sph':
 if geo_p=='ray':
  chi_nouv=snapshot_compl_per(r_nouv,r_v_0,config,res_gmsh,typ_sol)
elif config=='cylsph':
 if geo_p=='ray_sph':
  chi_nouv=snapshot_compl_per(r_nouv,r_c_0,config,res_gmsh,typ_sol)
 elif geo_p=='ray_cyl':
  chi_nouv=snapshot_compl_per(r_s_0,r_nouv,config,res_gmsh,typ_sol)

## Exploitation du champ ainsi obtenu
rho=r_nouv
r=r_nouv

plot(chi_nouv, linewidth=lw)
plt.title("Rho = 0,"+str(int(round(100*r_nouv,2))),fontsize=40)
if fig_todo=='aff':
 plt.show()
elif fig_todo=='save':
 plt.savefig("Figures3D/sol_rom"+str(int(round(100*r_nouv,2)))+"_sur"+str(Nsnap)+config+'_'+geo_p+"res"+str(res)+".png")
plt.close()
#sys.exit()
# Affichage des valeurs et erreurs de la solution périodique, quelle que soit la configuration
if config=='sph_un' or config=='cyl_un':
 #err_per_ind_01(chi_n,cen,r,npas_err)
 err_per_gr(cen_snap_ray,r_nouv,chi_nouv,npas_err,fig_todo)
elif config=='2sph':
 err_per_gr_compl(config,r_v_0,chi_nouv,npas_err,fig_todo)
elif config=='cylsph':
 if geo_p=='ray_sph':
  err_per_gr_compl(config,r_c_0,chi_nouv,npas_err,fig_todo)
 elif geo_p=='ray_cyl':
  err_per_gr_compl(config,r_nouv,chi_nouv,npas_err,fig_todo)
 #elif geo_p=='ray_linked':

# Tenseur de diffusion homogénéisé
## Intégrale de chi sur le domaine fluide
T_chi=np.zeros((3,3))
for k in range(0,3):
 for l in range(0,3):
  T_chi[k,l]=assemble(grad(chi_nouv)[k,l]*dx)
## Intégrale de l'identité sur le domaine fluide : voir ce qui précède avec lea porosité
print('Noeuds :',V_nouv.dim())
print('Porosité :',por)
#print('tenseur : ',T_chi)
## Calcul et affichage du tenseur Dhom
Dhom_kMEF=D_k*(D+T_chi.T)
#print(('Tenseur Dhom_k',Dhom_k))
print('Coefficient Dhom_k11 '+conf_mess+', '+geo_mess+' valeur '+str(rho)+ ' MEF :',Dhom_kMEF[0,0])

## Comparaison

err_rel=100*(Dhom_kMOR[0,0]-Dhom_kMEF[0,0])/Dhom_kMEF[0,0]
print('Erreur relative MEF-MOR :', err_rel , ' pourcent')

## On enregistre et imprime le temps d'éxécution de SE4

end=time.time()

print('se4 faite ',end-start,' secondes')
