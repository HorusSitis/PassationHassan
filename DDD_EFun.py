### Paquets à importer ###

from fenics import *
from dolfin import *
from mshr import *
import matplotlib.pyplot as plt
import numpy as np
from math import sqrt
import sys
import os

# Calcul parallèle

##import multiprocessing

# Performances

import time

# Stockage d'objets python

import marshal as ma
import shelve as sh

# Fonctions maison

from DDD_fun_obj import *

## Paquets spécifiques à la 3d ##

from mpl_toolkits.mplot3d.axes3d import get_test_data
# This import registers the 3D projection, but is otherwise unused.
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
import mpl_toolkits.mplot3d.art3d as art3d
from matplotlib.patches import Circle, PathPatch

##########################################################
### ------------ Code à lire : conditions ------------ ###
##########################################################

res_gmsh=50
config='2sph'#'cyl_un'#'sph_un'#'cylsph'#
un_snap_done=True
mesh_todo=True


# nom de l'appareil utilisé pour générer les données enregistrées
computer='MECALAC_29x8'#'T1700_35x8'#

# paramètres pour l'éxécution : affichage, tests de périodicité etc
fig_todo='save'
typ_msh='gms'
D_k=1.0
npas_err=20
typ_sol="bic_cyr"#"default"#seulement si res=10##

# Résolution EF : calcul parallèle ou séquentiel, prise en compte de la résolution
gen_snap='seq'#'seq_par'#

# répertoire pour les résultats

repertoire_parent="Res3D/"

########################################################################################################################################
## Calculs pour la section 2 de l'article-rapport-chapitre [...] . Un seul jeu de paramètres par configuration. Etapes chronométrées. ##
########################################################################################################################################

tol=1e-10

xinf=0.0
yinf=0.0
zinf=0.0
xsup=1.0
ysup=1.0
zsup=1.0

#determiner le domaine fixe pour interpoler la solution

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

### inclusions simples
if config=='sph_un':
 conf_mess='sphère unique'
elif config=='cyl_un':
 conf_mess='cylindre unique'
### inclusions composées
elif config=='2sph':
 conf_mess='deux sphères'
elif config=='cylsph':
 conf_mess='un cylindre et une sphère'

## Sphère unique
#if geo_p=='ray':
cen_snap_ray=[0.5,0.5,0.5]
## Cylindre unique
#if geo_p=='ray':
axe_snap_ray=[0.5,0.5]

r_s_0=0.35
r_c_0=0.25
r_v_0=0.2

res=res_gmsh

## inclusions simples

if config=='sph_un':
 mesh_prefix="cubesphere_periodique_triangle"
 mesh_name=mesh_prefix+"_"+str(int(round(100*r_s_0,2)))+"sur"+str(res)
elif config=='cyl_un':
 mesh_prefix="cubecylindre_periodique_triangle"
 mesh_name=mesh_prefix+"_"+str(int(round(100*r_c_0,2)))+"sur"+str(res)
## configurations complexes
elif config=='2sph':
 mesh_prefix='cube'+config+'_periodique_triangle'
 mesh_name=mesh_prefix+'_'+str(int(round(100*r_s_0,2)))+str(int(round(100*r_v_0,2)))+"sur"+str(res)
elif config=='cylsph':
 mesh_prefix='cube'+config+'_periodique_triangle'
 mesh_name=mesh_prefix+'_'+str(int(round(100*r_c_0,2)))+str(int(round(100*r_s_0,2)))+"sur"+str(res)

# ------------------------- Maillages chronométrés ------------------------- #

if mesh_todo:
 ## répertoire des maillages
 os.chdir(os.getcwd() + "/maillages_per/"+str(dimension)+"D")
 ## Génération d'un fichier .geo ? On commence avec un fichier unique et on modifie geo_p dans le code avant de sauvegarder sous le nom courant.
 print(mesh_name)
 ## Visualisation du fichier .geo
 print("gmsh "+mesh_name+".geo")
 if fig_todo=='aff':
  os.system("gmsh "+mesh_name+".geo")
 ## Conversion en .msh
 start=time()
 os.system("gmsh -"+str(dimension)+" "+mesh_name+".geo")
 end=time()
 tps_1=end-start
 print("temps de génération du maillage : "+str(tps_1)+" secondes")
 ## Affichage du maillage obtenu
 if fig_todo=='aff':
  os.system("gmsh "+mesh_name+".msh")
 if not un_snap_done:
  ## Conversion en .xml avec dolfin pour FEniCS : pas si on veut charger une solution EF calculée auparavant
  start=time()
  os.system("dolfin-convert "+mesh_name+".msh "+mesh_name+".xml")
  end=time()
  tps_2=end-start
  print("temps de conversion du maillage : "+str(tps_2)+" secondes")
  print("temps total d'éxécution : "+str(tps_1+tps_2)+" secondes")
 ## retour au répertoire PassationHassan
 os.chdir("../../")

# ------------------------- Snapshots, conditionnellement ------------------------- #

if not un_snap_done:
 # -------- Calcul des snapshots, sous forme vectorielle, avec des étiquettes -------- #
 ### Génération séquentielle des snapshots, pour des tests de la méthode des éléments finis ###
 #elif gen_snap=='seq':
 start=time()
 if config=='sph_un':
  chi_un=snapshot_sph_per(cen_snap_ray,r_s_0,res_gmsh,typ_sol)
 elif config=='cyl_un':
  chi_un=snapshot_cyl_per(axe_snap_ray,r_c_0,res_gmsh,typ_sol)
 elif config=='2sph':
  chi_un=snapshot_compl_per(r_s_0,r_v_0,config,res_gmsh,typ_sol)
 elif config=='cylsph':
  chi_un=snapshot_compl_per(r_s_0,r_c_0,config,res_gmsh,typ_sol)
 end=time()
 print('temps EF : ',end-start,' secondes')
 ### Génération parallèle pour chaque snapshot, pour de gros maillages ###
 #elif gen_snap=='seq_par':
 # chi_un=[]
 # vectorisation de la solution à sauvegarder
 chi_un_v=chi_un.vector().get_local()
 # Liste des snapshots : sauvegarde, on précise l'identité de la machine qui a effectué le calcul
 v_name='Chi_un_'+config+'_'+"sur"+str(res)+'_'+computer
 # sauvegarde de la liste des solutions indexées calculées avec la méthode des éléments finis
 with sh.open(repertoire_parent+v_name) as v_sto:
  v_sto["maliste"] = chi_un_v
 # Matrice des snapshots : plus tard, voir l'étape II
else :
 v_name='Chi_un_'+config+'_'+"sur"+str(res)+'_'+computer
 with sh.open(repertoire_parent+v_name) as v_loa:
  chi_un_v = v_loa["maliste"]

# --------------------------------------------------------------------------------- #
print('Taille de la solution EF :',len(chi_un_v))
# Exploitation des solution du problème aux éléments finis


mesh=Mesh("maillages_per/"+str(dimension)+"D/"+mesh_name+".xml")
V_un=VectorFunctionSpace(mesh, 'P', 2, constrained_domain=PeriodicBoundary())
print('Taille du maillage importé :',V_un.dim())
# On restitue la forme fonctionnelle de la solution
chi_un=Function(V_un)
chi_un.vector().set_local(chi_un_v)

if res==10:
 lw=0.27
elif res==20:
 lw=0.15
elif res==50:
 lw=0.01

# Représentation graphique
plot(chi_un, linewidth=lw)
plt.tight_layout(pad=0)
if fig_todo=='aff':
 plt.show()
else:
 plt.savefig("Figures3D/sol_un"+config+'_'+"res"+str(res)+".png")
plt.close()

# Tenseur de diffusion homogénéisé
## Intégrale de khi sur le domaine fluide
T_chi=np.zeros((3,3))
for k in range(0,3):
 for l in range(0,3):
  T_chi[k,l]=assemble(grad(chi_un)[k,l]*dx)
#print(T_chi)
## Intégrale de l'identité sur le domaine fluide
if config=='sph_un':
 por=1-4/3*pi*r_s_0**3
elif config=='cyl_un':
 por=1-pi*r_c_0**2
elif config=='2sph':
 por=1-4/3*pi*(r_s_0**3+r_v_0**3)
elif config=='cylsph' :
 por=1-4/3*pi*r_s_0**3-pi*r_c_0**2
D=por*np.eye(3)
## Calcul et affichage du tenseur Dhom
Dhom_k=D_k*(D+T_chi.T)
#print(('Tenseur Dhom_k',Dhom_k))
#print("Noeuds",V_un.dim())
#print("Porosité :",por)
#print('Coefficient Dhom_k11EF, snapshot unique '+", "+conf_mess+" :",Dhom_k[0,0])
integ=assemble(chi_un[1]*dx)
print('Valeur moyenne : ',integ)


## Anisotropie
mod_diag=min(abs(Dhom_k[0,0]),abs(Dhom_k[1,1]),abs(Dhom_k[2,2]))
mod_ndiag=0
for i in range(0,3):
 for j in range(0,i):
  if abs(Dhom_k[i,j])>mod_ndiag:
   mod_ndiag=abs(Dhom_k[i,j])
 for j in range(i+1,2):
  if abs(Dhom_k[i,j])>mod_ndiag:
   mod_ndiag=abs(Dhom_k[i,j])
print("Anisotropie I- : ",max(abs(Dhom_k[0,0]),abs(Dhom_k[1,1]),abs(Dhom_k[2,2]))/mod_diag-1)
print("Anisotropie II- : ",mod_ndiag/mod_diag)
#
