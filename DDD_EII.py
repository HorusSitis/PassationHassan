#####################################################################################################################################
######################################### Etape II : extrapolation des cliches, domaine_fixe ########################################
#####################################################################################################################################

# maillage du domaine fixe

mesh_dir="maillages_per/3D/"

## inclusions simples ou rayons lies
if dom_fixe=="am":
 mesh_f_name=mesh_dir+"cube_periodique_triangle"+"_"+dom_fixe+"_sur"+str(res_gmsh)+"_fixe.xml"
## inclusions multiples, unique rayon variable
elif dom_fixe=="solid":
 mesh_fixe_prefix=mesh_dir+"cube"+config+"_periodique_triangle_"
 if config=='2sph':
  mesh_f_name=mesh_fixe_prefix+"fixe"+str(int(round(100*r_v_0,2)))+"sur"+str(res_gmsh)+".xml"
 elif config=='cylsph':
  ## rayon du cylindre aux aretes ou de la sphere centrale fixes a 0.15 ##
  if geo_p=='ray_sph':
   mesh_f_name=mesh_fixe_prefix+str(int(round(100*r_c_0,2)))+"fixe"+"sur"+str(res_gmsh)+".xml"
  elif geo_p=='ray_cyl':
   if fixe_comp=='cyl_sph':
    mesh_f_name=mesh_fixe_prefix+"fixe"+str(int(round(100*r_s_0,2)))+"sur"+str(res_gmsh)+".xml"
   elif fixe_comp=='sph_un':
    mesh_f_name=mesh_dir+"cubesphere_periodique_triangle_"+str(int(round(100*r_s_0,2)))+"sur"+str(res_gmsh)+".xml"
elif dom_fixe=="ray_min":
 fixe_comp=True#utilisation du domaine fixe avec annulation du rayon du cylindre dans le fichier general
 if config=='cylsph':
  if geo_p=='ray_sph':
   mesh_f_name=mesh_dir+"cube"+config+"_periodique_triangle_"+str(int(round(100*r_c_0,2)))+str(int(round(100*r_min,2)))+"sur"+str(res_gmsh)+".xml"
  elif geo_p=='ray_cyl':
   mesh_f_name=mesh_dir+"cube"+config+"_periodique_triangle_"+str(int(round(100*r_min,2)))+str(int(round(100*r_s_0,2)))+"sur"+str(res_gmsh)+".xml"

print("Maillage fixe :",mesh_f_name)
mesh_fixe=Mesh(mesh_f_name)

# fonctions test du domaine fixe

V_fixe=VectorFunctionSpace(mesh_fixe,'P',2,constrained_domain=PeriodicBoundary())

### ------------ Etapes reproduites : dependances directes de Main3D ------------ ###

# Chargement de la liste des snapshots physiques

l_name='Lchi_'+str(N_snap)+'_'+config+'_'+geo_p+'_'+"sur"+str(res)+'_'+ordo+'_'+computer
print(l_name)
with sh.open(repertoire_parent+l_name) as l_loa:
 list_chi_v = l_loa["maliste"]

#sys.exit()
# Extrapolation au domaine Omega_fixe : inclusion spherique de rayon 0.0001, chi_prime defini sur ce domaine

def extra_snap(n):
 r=0.05*n
 # chargement du snapshot courant
 chi_n_v=list_chi_v[n-1]
 # mise sous forme d'une fonction EF
 if config=='sph_un':
  mesh_name=mesh_dir+"cubesphere_periodique_triangle_"+str(int(round(100*r,2)))+"sur"+str(res)+".xml"
 elif config=='cyl_un':
  mesh_name=mesh_dir+"cubecylindre_periodique_triangle_"+str(int(round(100*r,2)))+"sur"+str(res)+".xml"
 else:
  mesh_prefix=mesh_dir+"cube"+config+"_periodique_triangle_"
  if config=='2sph':
   r_v=0.15
   mesh_name=mesh_prefix+str(int(round(100*r,2)))+str(int(round(100*r_v,2)))+"sur"+str(res)+".xml"
  elif config=='cylsph':
   ## rayon du cylindre aux aretes ou de la sphere centrale fixes a 0.15 ##
   if geo_p=='ray_sph':
    r_s=r
    r_c=0.15
    mesh_name=mesh_prefix+str(int(round(100*r_c,2)))+str(int(round(100*r_s,2)))+"sur"+str(res)+".xml"
   elif geo_p=='ray_cyl':
    r_s=0.15
    r_c=r
    mesh_name=mesh_prefix+str(int(round(100*r_c,2)))+str(int(round(100*r_s,2)))+"sur"+str(res)+".xml"
 # quelle que soit la configuration
 #print(mesh_name)
 mesh=Mesh(mesh_name)
 V_n=VectorFunctionSpace(mesh, 'P', 2, constrained_domain=PeriodicBoundary())
 chi_n=Function(V_n)
 chi_n.vector().set_local(chi_n_v)
 # extrapolation du snapshot au domaine fixe
 chi_n.set_allow_extrapolation(True)
 chi_n_prime=interpolate(chi_n,V_fixe)
 # verification sur chi_n
 ###
 ## on range le snapshot dans une liste
 #list_snap.append(chi_n_prime)
 chi_n_prime_v=chi_n_prime.vector().get_local()
 return([n,chi_n_prime_v])



# Constitution de la matrice des snapshots

if not exsnap_done:
 pool=multiprocessing.Pool(processes=8)
 list_chi_n_prime_v=pool.map(extra_snap,(n for n in range(1,1+N_snap)))
 # Remplissage de la matrice
 nb_noeuds=V_fixe.dim()
 Usnap=np.zeros((nb_noeuds,N_snap))
 for n in range(1,1+N_snap):
  for i in range(0,N_snap):
   if list_chi_n_prime_v[i][0]==n:
    Usnap[:,n-1]=list_chi_n_prime_v[i][1]
 # Stochage de la matrice des snapshots
 u_name='Usnap_'+dom_fixe+'_'+str(N_snap)+'_'+config+'_'+geo_p+'_'+"res"+str(res)+'_'+ordo+'_'+computer
 #
 with sh.open(repertoire_parent+u_name) as u_sto:
  u_sto["maliste"] = Usnap
 #sys.exit()#------------------------------------------------------------
else:
 # Chargement de la matrice des snapshots
 u_name='Usnap_'+dom_fixe+'_'+str(N_snap)+'_'+config+'_'+geo_p+'_'+"res"+str(res)+'_'+ordo+'_'+computer
 with sh.open(repertoire_parent+u_name) as u_loa:
  Usnap = u_loa["maliste"]


#print(Usnap[0:5,0:5])
#sys.exit("attendons un peu pour les representations graphiques post-traitement")
# Representations graphiques

list_snap=[]
for n in range(1,1+N_snap):
 chi_prime=Function(V_fixe)
 chi_prime.vector().set_local(Usnap[:,n-1])
 # remplissage de la liste de fonctions
 list_snap.append(chi_prime)

#cen=cen_snap_ray
for n in range(1,1+N_snap):
 chi_prime_n=list_snap[n-1]
 # Affichage des valeurs de la solution interpolee
 plot(chi_prime_n, linewidth=lw)
 plt.title("Snapshot "+str(n),fontsize=30)
 if fig_todo=='aff':
  plt.show()
 else:
  plt.savefig("Figures3D/snap_interp"+dom_fixe+"_"+str(n)+"_sur"+str(N_snap)+config+'_'+geo_p+"res"+str(res)+".png")
 plt.close()
 # Affichage des valeurs et erreurs de la solution periodique, quelle que soit la configuration
 #err_per_ind_01(chi_prime_n,cen,r,npas_err)
 r=n*0.05
 ### Debuggage : cylindre et sphere ###
 if test_Dhom and (config=='cylsph' or config=='2sph'):
  if geo_p=='ray':
   r_cen=r
   r_per=r_v_0
   por=1-4/3*pi*(r_cen**3+r_per**3)
   por_prime=1-4/3*pi*r_per**3
   mesh_postfixe=str(int(round(100*r_cen,2)))+str(int(round(100*r_per,2)))+"sur"+str(res)
  elif geo_p=='ray_sph':
   r_s=r
   r_c=r_c_0
   por=1-4/3*pi*r_s**3-pi*r_c**2
   por_prime=1-4/3*pi*0.05**3-pi*r_c**2
   mesh_postfixe=str(int(round(100*r_c,2)))+str(int(round(100*r_s,2)))+"sur"+str(res)
  elif geo_p=='ray_cyl':
   r_s=r_s_0
   r_c=r
   por=1-4/3*pi*r_s**3-pi*r_c**2
   por_prime=1-4/3*pi*r_s**3-pi*0.05**2
   mesh_postfixe=str(int(round(100*r_c,2)))+str(int(round(100*r_s,2)))+"sur"+str(res)
  ##
  mesh_name="cube"+config+"_periodique_triangle_"+mesh_postfixe
  print("Verification : maillage courant",mesh_name)
  mesh=Mesh("maillages_per/3D/"+mesh_name+".xml")
  ##
  V_n=VectorFunctionSpace(mesh, 'P', 2, constrained_domain=PeriodicBoundary())
  chi_n=Function(V_n)
  chi_prime_n.set_allow_extrapolation(True)
  chi_n=interpolate(chi_prime_n,V_n)
  ##
  T_chi=np.zeros((3,3))
  for k in range(0,3):
   for l in range(0,3):
    T_chi[k,l]=assemble(grad(chi_n)[k,l]*dx)
  ##
  if geo_p=='ray':
   class DomPhysFluide(SubDomain):
    def inside(self, x, on_boundary):
     return True if (x[0]**2+x[1]**2+x[2]**2>=r_cen**2) else False
  elif geo_p=='ray_sph':
   class DomPhysFluide(SubDomain):
    def inside(self, x, on_boundary):
     return True if (x[0]**2+x[1]**2+x[2]**2>=r_s**2) else False
  elif geo_p=='ray_cyl':
   class DomPhysFluide(SubDomain):
    def inside(self, x, on_boundary):
     return True if (x[0]**2+x[2]**2>=r_c**2 and (1-x[0])**2+x[2]**2>=r_c**2 and x[0]**2+(1-x[2])**2>=r_c**2 and (1-x[0])**2+(1-x[2])**2>=r_c**2) else False
  dom_courant=DomPhysFluide()
  subdomains=MeshFunction('size_t',mesh_fixe,mesh_fixe.topology().dim())
  subdomains.set_all(1)
  dom_courant.mark(subdomains,12829)
  dxf=Measure("dx", domain=mesh_fixe, subdomain_data=subdomains)
  T_chi_restr_prime=np.zeros((3,3))
  for k in range(0,3):
   for l in range(0,3):
    T_chi_restr_prime[k,l]=assemble(grad(chi_prime_n)[k,l]*dxf(12829))
  ##
  T_chi_prime=np.zeros((3,3))
  for k in range(0,3):
   for l in range(0,3):
    T_chi_prime[k,l]=assemble(grad(chi_prime_n)[k,l]*dxf)
  #sys.exit()#------------------------------------------------------------------------------------------------------------------------
  ## les trois coefficients a comparer
  ###
  D=por*np.eye(3)
  integr_k_postprime=D_k*(D+T_chi.T)
  Dhom_k_postprime=integr_k_postprime*(1/por_prime)#por en dimensionnel
  ###
  D=por*np.eye(3)
  integr_k_restr_prime=D_k*(D+T_chi_restr_prime.T)
  Dhom_k_restr_prime=integr_k_restr_prime*(1/por_prime)
  ###
  D_prime=por_prime*np.eye(3)
  integr_k_prime=D_k*(D_prime+T_chi_prime.T)
  Dhom_k_prime=integr_k_prime*(1/por_prime)
  ##
 elif test_Dhom and config=='cyl_un':
  mesh_name="cubecylindre_periodique_triangle_"+str(int(round(100*r,2)))+"sur"+str(res)
  mesh=Mesh("maillages_per/3D/"+mesh_name+".xml")
  ##
  V_n=VectorFunctionSpace(mesh, 'P', 2, constrained_domain=PeriodicBoundary())
  chi_n=Function(V_n)
  chi_prime_n.set_allow_extrapolation(True)
  chi_n=interpolate(chi_prime_n,V_n)
  ##
  T_chi=np.zeros((3,3))
  for k in range(0,3):
   for l in range(0,3):
    T_chi[k,l]=assemble(grad(chi_n)[k,l]*dx)
  por=1-pi*r**2
  D=por*np.eye(3)
  integr_k_postprime=D_k*(D+T_chi.T)
  Dhom_k_postprime=integr_k_postprime*(1/por_prime)#por en dimensionnel
 ## ------------------------ Sphere unique, test ------------------------ ##
 elif test_Dhom and config=='sph_un':
  if geo_p=='ray':
   r_cen=r
   por=1-4/3*pi*r_cen**3
   por_prime=1
  ##
  mesh_name="cubesphere_periodique_triangle_"+str(int(round(100*r_cen,2)))+"sur"+str(res)
  print("Verification : maillage courant",mesh_name)
  mesh=Mesh("maillages_per/3D/"+mesh_name+".xml")
  ##
  V_n=VectorFunctionSpace(mesh, 'P', 2, constrained_domain=PeriodicBoundary())
  chi_n=Function(V_n)
  chi_prime_n.set_allow_extrapolation(True)
  chi_n=interpolate(chi_prime_n,V_n)
  ##
  T_chi=np.zeros((3,3))
  for k in range(0,3):
   for l in range(0,3):
    T_chi[k,l]=assemble(grad(chi_n)[k,l]*dx)
  ##
  class DomPhysFluide(SubDomain):
   def inside(self, x, on_boundary):
    return True if (x[0]**2+x[1]**2+x[2]**2>=r_cen**2) else False
  dom_courant=DomPhysFluide()
  subdomains=MeshFunction('size_t',mesh_fixe,mesh_fixe.topology().dim())
  subdomains.set_all(1)
  dom_courant.mark(subdomains,12829)
  dxf=Measure("dx", domain=mesh_fixe, subdomain_data=subdomains)
  T_chi_restr_prime=np.zeros((3,3))
  for k in range(0,3):
   for l in range(0,3):
    T_chi_restr_prime[k,l]=assemble(grad(chi_prime_n)[k,l]*dxf(12829))
  ##
  T_chi_prime=np.zeros((3,3))
  for k in range(0,3):
   for l in range(0,3):
    T_chi_prime[k,l]=assemble(grad(chi_prime_n)[k,l]*dxf)
  ## les trois coefficients a comparer
  ###
  D=por*np.eye(3)
  integr_k_postprime=D_k*(D+T_chi.T)
  Dhom_k_postprime=integr_k_postprime*(1/por_prime)#por en dimensionnel
  ###
  D=por*np.eye(3)
  integr_k_restr_prime=D_k*(D+T_chi_restr_prime.T)
  Dhom_k_restr_prime=integr_k_restr_prime*(1/por_prime)
  ###
  D_prime=por_prime*np.eye(3)
  integr_k_prime=D_k*(D_prime+T_chi_prime.T)
  Dhom_k_prime=integr_k_prime*(1/por_prime)
  ##
 if test_Dhom:
  print('Geometrie : '+conf_mess+', '+geo_mess+', '+str(int(round(100*r,2))))
  print('DUsnap fixe :',Dhom_k_prime[0,0],'porosite fixe',por_prime)
  print('DUsnap fixe restreint au domaine courant :',Dhom_k_restr_prime[0,0],'porosite',por)
  print('DUsnap physique :',Dhom_k_postprime[0,0],'porosite',por)
  print()
 ###
 # Dans ce cas, les faces du cube sont entieres
 ##cen=[0.5,0.5,0.5]
 ##err_per_gr(cen,r,chi_prime_n,npas_err,fig_todo)
 if config=='cylsph' and geo_p=='ray_sph':
  err_per_gr_compl(config,0.15,chi_prime_n,npas_err,fig_todo)
