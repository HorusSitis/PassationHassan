####################################################################################################################################
## Etape III : en utilisant la méthode des snapshots, calcul de la POD et des coefficients aléatoires, toujours dans domaine_fixe ##
####################################################################################################################################

#import PO23D as pod
#pod = reload(pod)
from PO23D import *
#
# Domaine d'interpolation et matrice des snapshots
#
c_x=0.5
c_y=0.5
c_z=0.5
#
#res=15
#
#
#Nsnap=8
#
#V_fixe=VectorFunctionSpace(mesh_fixe, "P", 2, constrained_domain=PeriodicBoundary())
nb_noeuds = V_fixe.dim()
#
Usnap=np.zeros((nb_noeuds,Nsnap))
#
pas=0.05
def inter_snap_ray(n):
 r=n*pas
 # Cliché sur le domaine avec inclusion
 u=snapshot_sph_per([c_x,c_y,c_z],r,res)
 # Extrapolation du cliché : khi prime
 u.set_allow_extrapolation(True)
 u_fixe=interpolate(u,V_fixe)
 # Forme vectorielle de la solution EF
 u_fixe_v=u_fixe.vector().get_local()
 return([n,u_fixe_v])
#
## Remplissage séquentiel de la matrice des snapshots
if rempUsnap=='seq':
 for n in range(1,1+Nsnap):
  u_fixe=inter_snap_ray(n)[1]
  # Remplissage de la matrice des snapshots
  Usnap[:,n-1]=u_fixe.vector().get_local()
## Remplissage parallèle de la matrice des snapshots
elif rempUsnap=='par8':
 pool=multiprocessing.Pool(processes=8)
 #
 Uf_par=pool.map(inter_snap_ray,(n for n in range(1,1+Nsnap)))
 #
 for n in range(1,1+Nsnap):
  for i in range(0,Nsnap):
   if Uf_par[i][0]==n:
    u_fixe_v=Uf_par[i][1]
    Usnap[:,n-1]=u_fixe_v
#
## matrice de corrélation
C=mat_corr_temp(V_fixe,Nsnap,Usnap)
#
## Calcul des coefficients aléatoires et la base POD
vp_A_phi=mat_a_mat_phi(Nsnap,Usnap,C)
#
val_propres=vp_A_phi[0]
Aleat=vp_A_phi[1]
phi=vp_A_phi[2]
print(val_propres)
#res, res_fixe=20 : énergie [71%, 24%, 5%, 0.37%, 0.058%, 0%, 0%, 0%] 
