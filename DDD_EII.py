#####################################################################################################################################
######################################### Etape II : extrapolation des clichés, domaine_fixe ########################################
#####################################################################################################################################

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

# Chargement de la liste des snapshots physiques

l_name='Lchi_'+str(Nsnap)+'_'+config+'_'+geo_p+'_'+ordo+'_'+computer

with sh.open(repertoire_parent+l_name) as l_loa:
    list_chi_v = l_loa["maliste"]

# Extrapolation au domaine Omega_fixe : inclusion sphérique de rayon 0.0001, chi_prime défini sur ce domaine

mesh_fixe=Mesh("maillages_per/3D/cubesphere_periodique_triangle_0001fixe.xml")

V_fixe=V=VectorFunctionSpace(mesh_fixe, 'P', 2, constrained_domain=PeriodicBoundary())

list_snap=[]

for n in range(1,1+Nsnap):
 # chargement du snapshot courant
 chi_n_v=list_chi_v[n-1]
 r=0.05*n
 if typ_msh=='gms':
  mesh=Mesh("maillages_per/3D/cubesphere_periodique_triangle_"+str(int(round(100*r,2)))+".xml")
 #else:mesh=creer_maill_sph(cen,r,res)
 V_n=VectorFunctionSpace(mesh, 'P', 2, constrained_domain=PeriodicBoundary())
 chi_n=Function(V_n)
 chi_n.vector().set_local(chi_n_v)
 # extrapolation du snapshot au domaine fixe
 chi_n.set_allow_extrapolation(True)
 chi_n_prime=interpolate(chi_n,V_fixe)
 # on range le snapshot dans une liste
 list_snap.append(chi_n_prime)
# Stockage des snapshots virtuels obtenus

#l_name='Lchi_extra_'+str(Nsnap)+'_'+ordo+'_'+computer

#with sh.open(repertoire_parent+l_name) as l_sto:
#    l_sto["maliste"] = list_snap

# Constitution de la matrice des snapshots


nb_noeuds=V_fixe.dim()

Usnap=np.zeros((nb_noeuds,Nsnap))

for n in range(0,Nsnap):
 chi_n_prime=list_snap[n]
 Usnap[:,n]=chi_n_prime.vector().get_local()

# Stochage de la matrice des snapshots

u_name='Usnap_'+str(Nsnap)+'_'+config+'_'+geo_p+'_'+ordo+'_'+computer

with sh.open(repertoire_parent+u_name) as u_sto:
    u_sto["maliste"] = Usnap

# Représentations graphiques

cen=cen_snap_ray
for n in range(1,1+Nsnap):
 chi_prime_n=list_snap[n-1]
 # Affichage des valeurs de la solution interpolée
 plot(chi_prime_n)
 if fig_todo=='aff':
  plt.show()
 else:
  plt.savefig("Figures3D/snap_interp0001_"+str(n)+"_sur"+str(Nsnap)+config+geo_p+".png")
 plt.close()
 # Affichage des valeurs et erreurs de la solution périodique, quelle que soit la configuration
 #err_per_ind_01(chi_prime_n,cen,r,npas_err)
 err_per_gr(cen,r,chi_prime_n,npas_err,fig_todo)
