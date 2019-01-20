#################################################################################################
## Etape IV : Prédictions. Choisir les paramètres du problème à résoudre par le modèle réduit. ##
#################################################################################################

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

#r_nouv=0.22

if typ_msh=='gms':
 mesh=Mesh("maillages_per/3D/cubesphere_periodique_triangle_"+str(int(round(100*r_nouv,2)))+".xml")
#else:mesh=creer_maill_sph(cen,r_nouv,res)

V_nouv=VectorFunctionSpace(mesh_fixe, "P", 2, constrained_domain=PeriodicBoundary())

# --------------------- SE1 : projection de la base POD sur le nouveau domaine --------------------- #

mesh_fixe=Mesh("maillages_per/3D/cubesphere_periodique_0001fixe.xml")
V_fixe=V=VectorFunctionSpace(mesh_fixe, 'P', 2, constrained_domain=PeriodicBoundary())

## Chargement de la marice des snapshots

u_name='Usnap_'+str(Nsnap)+'_'+config+'_'+geo_p+'_'+ordo+'_'+computer

with sh.open(repertoire_parent+l_name) as u_loa:
    Usnap = u_loa["maliste"]

## Choix de la base pour le modèle réduit




## Snapshots sous forme fonctionnelle, extrapolation au domaine V_nouv

nb_noeuds_nouv=V_nouv.dim()
Phi_nouv=np.zeros((nb_noeuds_nouv,nb_modes))






# --------------------- SE2 : résolution du modèle réduit --------------------- #
