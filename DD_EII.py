#####################################################################################################################################
######################################### Etape II : extrapolation des clichés, domaine_fixe ########################################
#####################################################################################################################################

### ------------ Reproduire éventuellement pour des étapes ultérieures. Laisser seulement dans DD_fun_obj ? ------------ ###

import time

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
# Chargement de la liste des snapshots physiques

l_name='Lchi_'+str(Nsnap)+'_'+config+'_'+geo_p+'_deg'+str(VFS_degree)+'_'+ordo+'_'+computer

with sh.open(repertoire_parent+l_name) as l_loa:
 list_chi_v = l_loa["maliste"]

# Extrapolation au domaine Omega_fixe :
if dom_fixe=="am":
 mesh_fixe=Mesh("maillages_per/2D/maillage_fixe2D_am.xml")
elif config=='compl':
 mesh_fixe=Mesh("maillages_per/2D/maillage_trous2D_"+geo_p+"_fixe.xml")
elif dom_fixe=="ray_min":
 if config=='cer_un':
  mesh_fixe=Mesh('maillages_per/2D/maillage_trou2D_5.xml')

V_fixe=VectorFunctionSpace(mesh_fixe, 'P', VFS_degree, constrained_domain=PeriodicBoundary())
#sys.exit()
# Extrapolation des solutions du problème physique

def extra_snap(n):
 if geo_p=='hor':
  r=0.01+0.04*n
 else:
  r=0.05*n
 # chargement du snapshot courant
 chi_n_v=list_chi_v[n-1]
 # mise sous forme d'une fonction EF
 if config!='compl':
  if cen_snap_ray==[0.,0.]:
   mention="_som"
  else:
   mention=""
  if typ_msh=='gms':
   mesh_name="maillages_per/2D/maillage_trou2D"+mention+"_"+str(int(round(100*r,2)))+".xml"
   print(mesh_name)
   mesh=Mesh(mesh_name)
 else:
  mesh_name="maillages_per/2D/maillage_trous2D_"+geo_p+"_"+str(int(round(100*r,2)))+".xml"
  mesh=Mesh(mesh_name)
 V_n=VectorFunctionSpace(mesh, 'P', VFS_degree, constrained_domain=PeriodicBoundary())
 chi_n=Function(V_n)
 chi_n.vector().set_local(chi_n_v)
 # extrapolation du snapshot au domaine fixe
 chi_n.set_allow_extrapolation(True)
 chi_n_prime=interpolate(chi_n,V_fixe)
 ## on range le snapshot dans une liste
 #list_snap.append(chi_n_prime)
 chi_n_prime_v=chi_n_prime.vector().get_local()
 return([n,chi_n_prime_v])


#cc=extra_snap(1)
#nb_noeuds=V_fixe.dim()
#snap=np.zeros((nb_noeuds,1))#Nsnap))
#print(len(cc[1]),nb_noeuds)

#sys.exit()

# Constitution de la matrice des snapshots

if not exsnap_done:
 pool=multiprocessing.Pool(processes=8)
 list_chi_n_prime_v=pool.map(extra_snap,(n for n in range(1,1+Nsnap)))
 # Remplissage de la matrice
 nb_noeuds=V_fixe.dim()
 Usnap=np.zeros((nb_noeuds,Nsnap))
 for n in range(1,1+Nsnap):
  for i in range(0,Nsnap):
   if list_chi_n_prime_v[i][0]==n:
    Usnap[:,n-1]=list_chi_n_prime_v[i][1]
 # Stochage de la matrice des snapshots
 u_name='Usnap_'+dom_fixe+str(Nsnap)+'_'+config+'_'+geo_p+'_deg'+str(VFS_degree)+'_'+ordo+'_'+computer
 #
 with sh.open(repertoire_parent+u_name) as u_sto:
  u_sto["maliste"] = Usnap
else:
 # Chargement de la matrice des snapshots
 u_name='Usnap_'+dom_fixe+str(Nsnap)+'_'+config+'_'+geo_p+'_deg'+str(VFS_degree)+'_'+ordo+'_'+computer
 print(repertoire_parent+u_name)
 with sh.open(repertoire_parent+u_name) as u_loa:
  Usnap = u_loa["maliste"]



# Représentations graphiques

list_snap=[]
for n in range(1,1+Nsnap):
 chi_prime=Function(V_fixe)
 chi_prime.vector().set_local(Usnap[:,n-1])
 # remplissage de la liste de fonctions
 list_snap.append(chi_prime)

cen=cen_snap_ray
for n in range(1,1+Nsnap):
 if geo_p=='hor':
  r=0.01+0.04*n
 else:
  r=0.05*n
 chi_prime_n=list_snap[n-1]
 if fig_todo=='aff' or fig_todo=='save':
  # Affichage des valeurs de la solution interpolée
  plot(chi_prime_n)
  plt.title("Snapshot "+str(n),fontsize=40)
  if fig_todo=='aff':
   plt.show()
  else:
   plt.savefig("Figures2D/snap_"+str(n)+"_sur"+str(Nsnap)+config+'_'+geo_p+".png")
  plt.close()
 else:
  print('pffrrh !')
 # Affichage des valeurs et erreurs de la solution périodique, quelle que soit la configuration
 #err_per_ind_01(chi_prime_n,cen,r,npas_err)
 #err_per_gr(cen,r,chi_prime_n,npas_err,fig_todo)
 # Affichage des composantes scalaires : interpolée
 if config=='cer_un' and geo_p=='ray' and fig_todo!='':
  fig_chi(cen_snap_ray,r,chi_prime_n,fig_todo)

if not test_Dhom:
 sys.exit('fin de l étape II, sans tests d intégration sur des sousdomaines')#-------------------------------------------------------------------------------------------------------------------------------------------

lg_crow=-1
crow=10**(lg_crow)
Nrefine=4
n_mp_refi=8

def f_testDhom(n):
 if geo_p=='hor':
  r=0.01+0.04*n
 else:
  r=0.05*n
 chi_prime_n=list_snap[n-1]
 if geo_p=='ray':
  por=1-pi*r**2
 por_prime=1
 ## Création du maillage fixe local
 if dom_fixe=="am":
  mesh_fixe=Mesh("maillages_per/2D/maillage_fixe2D_am.xml")
 elif config=='compl':
  mesh_fixe=Mesh("maillages_per/2D/maillage_trous2D_"+geo_p+"_fixe.xml")
 elif dom_fixe=="ray_min":
  if config=='cer_un':
   mesh_fixe=Mesh('maillages_per/2D/maillage_trou2D_5.xml')
 V_fixe=VectorFunctionSpace(mesh_fixe, 'P', VFS_degree, constrained_domain=PeriodicBoundary())
 ## Interpolation post-snapshot sur le maillage physique
 mesh_name="maillages_per/2D/maillage_trou2D"+mention+"_"+str(int(round(100*r,2)))+".xml"
 mesh=Mesh(mesh_name)
 V_n=VectorFunctionSpace(mesh,'P',VFS_degree,constrained_domain=PeriodicBoundary())
 chi_postprime_n=Function(V_n)
 chi_prime_n.set_allow_extrapolation(True)
 chi_postprime_n=interpolate(chi_prime_n,V_n)
 T_chi_postprime=np.zeros((2,2))
 for k in range(0,2):
  for l in range(0,2):
   T_chi_postprime[k,l]=assemble(grad(chi_postprime_n)[k,l]*dx)
 ## Intégration sur le domaine fluide avec le maillage fixe
 # Raffinement du maillage
 start=time.time()
 for i in range(Nrefine):
  print('raffinement',i)
  markers = MeshFunction("bool", mesh_fixe, mesh_fixe.topology().dim())
  markers.set_all(False)
  for c in cells(mesh_fixe):
   for f in facets(c):
    if ((f.midpoint()[0]-cen_snap_ray[0])**2+(f.midpoint()[1]-cen_snap_ray[1])**2<=(r*(1+crow))**2) and ((f.midpoint()[0]-cen_snap_ray[0])**2+(f.midpoint()[1]-cen_snap_ray[1])**2>=(r*(1-crow))**2):
     markers[c]=True
   for v in vertices(c):
    if (v.point().x()-cen_snap_ray[0])**2+(v.point().y()-cen_snap_ray[1])**2>=(r*(1-crow))**2 and (v.point().x()-cen_snap_ray[0])**2+(v.point().y()-cen_snap_ray[1])**2<=(r*(1+crow))**2:
     markers[c]=True
  mesh_fixe=refine(mesh_fixe, markers, redistribute=True)
 end=time.time()
 #print('Raffinemenent du maillage :',end-start,'secondes')
 tps_refi=end-start 
 # Interpolation sur le maillage raffiné
 if Nrefine>0:
  print('interpolation rayon',str(int(round(100*r,2))))
  V_fixe=VectorFunctionSpace(mesh_fixe, 'P', VFS_degree, constrained_domain=PeriodicBoundary())
  start=time.time()
  chi_prime_n.set_allow_extrapolation(True)
  chi_prime_n=interpolate(chi_prime_n,V_fixe)
  end=time.time()
  tps_interp=end-start
  #print('Interpolation sur le maillage raffiné :',end-start,'secondes')
 plot(mesh_fixe)
 plt.title('Maillage raffiné '+str(Nrefine)+' fois rayon '+str(int(round(100*r,2)))+'x10e-2')
 plt.show()
 plt.close()
 # Création du domaine d'intégration sur le maillage fixe
 class DomPhysFluide(SubDomain):
  def inside(self, x, on_boundary):
   return True if (x[0]**2+x[1]**2>=r**2) else False
 dom_courant=DomPhysFluide()
 subdomains=MeshFunction('size_t',mesh_fixe,mesh_fixe.topology().dim())
 subdomains.set_all(1)
 dom_courant.mark(subdomains,12829)
 dxf=Measure("dx", domain=mesh_fixe, subdomain_data=subdomains)
 T_chi_restr_prime=np.zeros((2,2))
 co_T_chi_restr_prime=np.zeros((2,2))
 for k in range(0,2):
  for l in range(0,2):
   T_chi_restr_prime[k,l]=assemble(grad(chi_prime_n)[k,l]*dxf(12829))
   co_T_chi_restr_prime[k,l]=assemble(grad(chi_prime_n)[k,l]*dxf(1))
 return([tps_refi,tps_interp,T_chi_restr_prime,co_T_chi_restr_prime,T_chi_postprime])

pool=multiprocessing.Pool(processes=n_mp_refi)
list_refi_interp=pool.map(f_testDhom,(n for n in range(1,1+Nsnap)))

nom_fichier='Res2D/'+'testDhom'+config+geo_p+'_'+dom_fixe+'_crown10'+str(lg_crow)+'_Nrefi'+str(Nrefine)+'.txt'
registre=open(nom_fichier,'w')

for n in range(1,1+Nsnap):
 if geo_p=='hor':
  r=0.01+0.04*n
 else:
  r=0.05*n
 if geo_p=='ray':
  por=1-pi*r**2
 por_prime=1
 ## Interpolation post-snapshot sur le maillage physique
 mesh_name="maillages_per/2D/maillage_trou2D"+mention+"_"+str(int(round(100*r,2)))+".xml"
 mesh=Mesh(mesh_name)
 V_n=VectorFunctionSpace(mesh,'P',VFS_degree,constrained_domain=PeriodicBoundary())
 chi_postprime_n=Function(V_n)
 chi_prime_n.set_allow_extrapolation(True)
 chi_postprime_n=interpolate(chi_prime_n,V_n)
 T_chi_postprime=np.zeros((2,2))
 for k in range(0,2):
  for l in range(0,2):
   T_chi_postprime[k,l]=assemble(grad(chi_postprime_n)[k,l]*dx)
 ## Chargement des résultats d'intégration des tenseurs d'homogénéisation sur les différents maillages
 [tps_refi,tps_interp,T_chi_restr_prime,co_T_chi_restr_prime,T_chi_postprime]=list_refi_interp[n-1]
 ## Intégration sur le domaine virtuel : carré unité éventuellement privé d'un point
 T_chi_prime=np.zeros((2,2))
 for k in range(0,2):
  for l in range(0,2):
   T_chi_prime[k,l]=assemble(grad(chi_prime_n)[k,l]*dx)
 ## les trois coefficients à comparer
 ###
 D=por*np.eye(2)
 integr_k_postprime=D_k*(D+T_chi_postprime.T)
 Dhom_k_postprime=integr_k_postprime*(1/por_prime)#por en dimensionnel
 ###
 integr_k_restr_prime=D_k*(D+T_chi_restr_prime.T)
 Dhom_k_restr_prime=integr_k_restr_prime*(1/por_prime)
 ###
 D_moy=D_k*(D+0.5*(T_chi_restr_prime.T-co_T_chi_restr_prime.T))*(1/por_prime)
 ###
 D_prime=por_prime*np.eye(2)
 integr_k_prime=D_k*(D_prime+T_chi_prime.T)
 Dhom_k_prime=integr_k_prime*(1/por_prime)
 ##
 print('##################################################################')
 print('Géométrie : '+conf_mess+', '+geo_mess+', '+str(int(round(100*r,2))))
 print('Domaine fixe : '+dom_fixe)
 print('------------------------------------------------------------------')
 print('Couronne :',crow)
 print('Nombre de tours de raffinement :',Nrefine)
 print('Raffinemenent du maillage :',tps_refi,'secondes')
 print('Interpolation sur le maillage raffiné :',tps_interp,'secondes')
 print('------------------------------------------------------------------')
 #print('DUsnap fixe :',Dhom_k_prime[0,0],'porosité fixe',por_prime)
 #print('Tenseur virtuel :',T_chi_prime)
 #print('Tenseur sur le domaine fluide :',T_chi_restr_prime)
 #print('Complémentaire :',co_T_chi_restr_prime)
 print()
 print('DUsnap fixe restreint au domaine courant :',Dhom_k_restr_prime[0,0])#,'porosité',por)
 print('DUsnap physique :',Dhom_k_postprime[0,0])#,'porosité',por)
 print('DUsnap domaine fixe moyenné :',D_moy[0,0])
 print('------------------------------------------------------------------')
 print('Erreur relative :',100*(Dhom_k_restr_prime[0,0]-Dhom_k_postprime[0,0])/Dhom_k_postprime[0,0],'pourcent')
 print('##################################################################')
 #
 registre.write('##################################################################'+'\n')
 registre.write('Rho : '+str(int(round(100*r,2)))+'\n')
 registre.write('Domaine fixe : '+dom_fixe+'\n')
 registre.write('------------------------------------------------------------------'+'\n')
 registre.write('Couronne : '+str(crow)+'\n')
 registre.write('Nombre de tours de raffinement : '+str(Nrefine)+'\n')
 registre.write('Raffinemenent du maillage : '+str(tps_refi)+' secondes'+'\n')
 registre.write('Interpolation sur le maillage raffine : '+str(tps_interp)+' secondes'+'\n')
 registre.write('------------------------------------------------------------------'+'\n')
 #registre.write('DUsnap fixe : '+str(Dhom_k_prime[0,0],'porosité fixe'+str(por_prime)+'\n')
 registre.write('Tenseur virtuel : '+str(T_chi_prime)+'\n')
 registre.write('DUsnap fixe restreint au domaine courant : '+str(Dhom_k_restr_prime[0,0])+'\n')
 registre.write('DUsnap physique : '+str(Dhom_k_postprime[0,0])+'\n')
 registre.write('DUsnap domaine fixe moyenne : '+str(D_moy[0,0])+'\n')
 registre.write('------------------------------------------------------------------'+'\n')
 registre.write('Erreur relative : '+str(100*(Dhom_k_restr_prime[0,0]-Dhom_k_postprime[0,0])/Dhom_k_postprime[0,0])+' pourcent'+'\n')
 registre.write('##################################################################'+'\n')


registre.close()
