### Répertoire courant ###

cd /home/amorea12/Documents/T_LaSIE/PassationHassan

### Paquets à importer ###

##############################################################################################################################
############################### Calculs avec la POD, modèles réduits ; merci à Hassan GHRAIEB. ###############################
##############################################################################################################################

#Paquets spécifiques à POD-MOR

from fenics import *
from dolfin import *
from mshr import *
import matplotlib.pyplot as plt
import numpy as np
from math import sqrt
import sys

import multiprocessing

#from DDD_fun_obj import *
import DDD_fun_obj as F3d
from importlib import reload
F3d=reload(F3d)

##############################################################################################################################
##############################################################################################################################
##############################################################################################################################

## Etape I : réalisation des clichés, avec la méthode des éléments finis. Stockage dans snap2D/ ##
# Utilise fenics, éventuellement augmenté par dolfin #

tol=1e-10

xinf=0.0
yinf=0.0
zinf=0.0
xsup=1.0
ysup=1.0
zsup=1.0

#determiner le domaine fixe pour interpoler la solution

class PeriodicBoundary(SubDomain):
 # Left boundary is "target domain" G
 def inside(self, x, on_boundary):
  return on_boundary and (near(x[0],xinf,tol) or near(x[1],yinf,tol) or near(x[2],zinf,tol))
 # Map right boundary (H) to left boundary (G)
 def map(self, x, y):
  if (near(x[0],xsup,tol)):
   y[0] = x[0] - 1.0
   y[1] = x[1]
   y[2] = x[2]        
  elif (near(x[0],xsup,tol)):
   y[0] = x[0]
   y[1] = x[1] - 1.0
   y[2] = x[2]
  else:
   y[0] = x[0]
   y[1] = x[1]
   y[2] = x[2] - 1.0

### Maillage sur le domaine Omega_fixe : aucune inclusion, porosité égale à 1 ###

res_fixe=15

domaine_fixe=Box(Point(xinf,yinf,zinf),Point(xsup,ysup,zsup))
mesh_fixe=generate_mesh(domaine_fixe,res_fixe)
V_fixe=VectorFunctionSpace(mesh_fixe, "P", 2, constrained_domain=PeriodicBoundary())

#représentation graphique du maillage
plot(mesh_fixe)
plt.show()
plt.close()

### Etape 1 : snapshots. Tests. ###

c_x=0.5
c_y=0.5
c_z=0.5

r=0.3

res=10

# Maillage raffiné unique #

mesh_c_r=F3d.creer_maill_sph([c_x,c_y,c_z],r,res)

plot(mesh_c_r)
plt.show()
plt.close()

# Champ khi unique #

res=20

u=F3d.snapshot_sph_per([c_x,c_y,c_z],r,res)

plot(u)
plt.show()
plt.close()

## Erreur de périodicité ##

# Boucles pour une famille de snapshots : variation du rayon ou du centre #


## Etape II : extrapolation des clichés, domaine_fixe ##

# Attention aux codes de Hassan : erreurs ... visibles sur les figures du rapport, s'il s'agit bien des snapshots extrapolés

c_x=0.5
c_y=0.5
c_z=0.5

res=15

for i in range(1,2):
 r=i*0.2#05
 u=F3d.snapshot_sph_per([c_x,c_y,c_z],r,res)
 ## chargement du snapshot pour l'indice courant
 # Extrapolation au domaine Omega_fixe : aucune inclusion, khi défini sur [0,1]times[0,1]
 u.set_allow_extrapolation(True)
 u_fixe=interpolate(u,V_fixe)##rapide
 u_fixe = project(u, V_fixe)##lent
 plot(u_fixe)
 plt.show()
 plt.close()

## Etape III : en utilisant la méthode des snapshots, calcul de la POD et des coefficients aléatoires, toujours dans domaine_fixe ##

import Antoine_POD as pod
pod = reload(pod)
#from Antoine_POD import *

# Domaine d'interpolation et matrice des snapshots

c_x=0.5
c_y=0.5
c_z=0.5

res=15

Nsnap=8

V_fixe=VectorFunctionSpace(mesh_fixe, "P", 2, constrained_domain=PeriodicBoundary())
nb_noeuds = V_fixe.dim()

Usnap=np.zeros((nb_noeuds,Nsnap))

pas=0.05
def inter_snap_ray(n):
 r=n*pas
 # Cliché sur le domaine avec inclusion
 u=F3d.snapshot_sph_per([c_x,c_y,c_z],r,res)
 # Extrapolation du cliché : khi prime
 u.set_allow_extrapolation(True)
 u_fixe=interpolate(u,V_fixe)
 # Forme vectorielle de la solution EF
 u_fixe_v=u_fixe.vector().get_local()
 return([n,u_fixe_v])

# Remplissage séquentiel de la matrice des snapshots
for n in range(1,1+Nsnap):
 u_fixe=inter_snap_ray(n)[1]
 # Remplissage de la matrice des snapshots
 Usnap[:,n-1]=u_fixe.vector().get_local()

# Remplissage paralèle de la matrice des snapshots
pool=multiprocessing.Pool(processes=8)

Uf_par=pool.map(inter_snap_ray,(n for n in range(1,1+Nsnap)))

for n in range(1,1+Nsnap):
 for i in range(0,Nsnap):
  if Uf_par[i][0]==n:
   u_fixe_v=Uf_par[i][1]
   Usnap[:,n-1]=u_fixe_v#.vector().get_local()

# matrice de corrélation
C=pod.mat_corr_temp(V_fixe,Nsnap,Usnap)

# Calcul des coefficients aléatoires et la base POD
vp_A_phi=pod.mat_a_mat_phi(Nsnap,Usnap,C)

val_propres=vp_A_phi[0]
Aleat=vp_A_phi[1]
phi=vp_A_phi[2]



