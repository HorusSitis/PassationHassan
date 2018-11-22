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
#import numpy as np
from math import sqrt
import sys

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

### Etape 1 : snapshots. Tests. ###

c_x=0.15
c_y=0.15
c_z=0.15

r=0.4

res=10

# Maillage raffiné unique #

mesh_c_r=F3d.creer_maill_sph([c_x,c_y,c_z],r,res)

plot(mesh_c_r)
plt.show()
plt.close()

# Champ khi unique #

u=F3d.snapshot_sph_per([c_x,c_y,c_z],r,res)

plot(u)
plt.show()
plt.close()

## Interpolation du champ ##

## Erreur de périodicité ##







# Boucles pour une famille de snapshots : variation du rayon ou du centre #




