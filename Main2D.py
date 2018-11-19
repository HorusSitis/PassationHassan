### Répertoire courant ###

cd /home/amorea12/Documents/T_LaSIE/PassationHassan

### Paquets à importer ###

import numpy as np
import random as rd

import pylab as pl
from pylab import *
#import matplotlib.pylab as pl

import os

#threading

import threading
import time

#multiprocessing


import multiprocessing
from joblib import Parallel, delayed

#stockage d'objets python

import marshal as ma

import shelve as sh

### Algorithme RSAA arrêté ###

## Codes à importer

from importlib import reload

###############################################################################################################
############################### Génération de structures périodiques aléatoires ###############################
###############################################################################################################

import RSAA_2d_ray as R2d

R2d = reload(R2d)

## Génération d'inclusions dans un réseau carré, stockage ##

cote=70#500 ; 1000
Tps=10000#1000000

cote=150
Tps=20000

cote=500
Tps=30000

cote=1000
Tps=1000000

# géométrie euclidienne

geo=[R2d.eucl2D,R2d.vol2D]

# lois normales pour la taille des inclusions ; paramètres #

Lois=[R2d.g_norm,R2d.g_norm]
Par=[(3,1),(1.5,5)]#[(5,0.5),(8,0.1)] pour 1000

# fractions volumiques pour les différentes phases ; marges

frac=[0.15,0.12]
delta_inc=2

L=R2d.RSAA_ph_dist2D(geo,delta_inc,Lois,Par,frac,[cote,cote],Tps,'per')

l_name='L_'+str(cote)+'_'+str(Tps)

rep_inc='IncAlea2D'

# stockage : listes de centres-rayons-phases et fractions volumiques #

#ma.dump(LLL, open("L_20", 'wb')) ## Sauvegarde de la liste 
#L_load = ma.load(open("L_20", "rb")) ## Rechargement de la liste

# sauvegarde de la variable maliste sous le nom "maliste" dans le fichier 'L_20_10000'
with sh.open(rep_inc+'/'+l_name) as l_sto:
    l_sto["maliste"] = L
 
# chargement de la variable maliste, on recharge les variables de dimension spatiale et du nombre de tours de boucle RSAA

#cote=
#Tps=
l_name='L_'+str(cote)+'_'+str(Tps)

with sh.open(rep_inc+'/'+l_name) as l_loa:
    L_loa = l_loa["maliste"]

## création d'une matrice contenant les inclustions : un coefficient pour un pixel ##

A_remp=R2d.Vremp2D(L_loa,R2d.eucl2D,[cote,cote],'per')

#pour avoir une fonction top-level à paralléliser

L=L_loa[0]
M=np.array(L)
n_phi=max(M[:,2])

dim=[cote,cote]
C_per='per'
dist=R2d.eucl2D

def parc_liste_k(k):
 return(R2d.parc_liste(k,L,n_phi,dist,dim,C_per))

pool=multiprocessing.Pool(processes=8)

B_par=pool.map(parc_liste_k,(k for k in range(dim[0]*dim[1])))

A_remp=np.array(B_par)
A_remp=np.reshape(A_remp,(dim[0],dim[1]))

##500 ; tps 30 000 ; 4600 inclusions : environ 1h sur MECALAC
##1000 ; tps 1 000 000 ; 2800 inclusions : environ 2 h 40 min pour l'ordinateur fixe 8 x 3,50GHz

#avec du calcul parallèle : imports avec joblib

#A_remp=R2d.ParVremp2D(L_loa,R2d.eucl2D,[cote,cote],'per')

##70 ; tps 10 000 ; ... inclusions : 20 secondes environ
##150 ; tps 20 000 ; 1000 inclusions : 2 minutes 20 secondes
##1000 ; tps 1 000 000 ; 2800 inclusions : 15h ++

# stockage de la matrice A #

a_name='VA_'+str(cote)+'_'+str(Tps)
#répertoire ?

rep_inc='IncAlea2D'

with sh.open(rep_inc+'/'+a_name) as a_sto:
    a_sto["maliste"] = A_remp

#cote=
#Tps=
a_name='VA_'+str(cote)+'_'+str(Tps)

with sh.open(rep_inc+'/'+a_name) as a_loa:
    A_loa = a_loa["maliste"]

## représentation graphique ##

# Couleurs : bleu pour le fluide, gris pour une première inclusion, jaune pour une deuxième
##avec pylab
cmap=matplotlib.colors.ListedColormap(['blue','grey','yellow'])
bounds=[0,1,2,3]
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
#

pl.imshow(A_remp,interpolation='none',cmap=cmap,norm=norm)
pl.axis('off')

figname=str(cote)+'f'+str(cote)+'Tps'+str(Tps)+'.png'
rep='Figures2D'
save_name=rep+'/'+figname

savefig(save_name)

pl.show()


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

#Import de codes de Hassan, pas de modification a priori. Reloads ?

cd 2DHassan/
#import calc_POD as PoD
#import Lecture_ecriture_homog as LEct
from calc_POD import *
from Lecture_ecriture_homog import *
#POD= reload(POD)
cd ../

#Imports divers : fonctions, objets etc

from DD_fun_obj import *#attention pas de chiffre au début du nom d'un paquet
r=0.0
## problème à l'import : r prend la valeur index(6)


import DD_fun_obj as fun_obj

fun_obj=reload(fun_obj)

### Codes éxécutés : cas d'une inclusion circulaire, le rayon du disque central est le paramètre pour la POD ###

## Etape I : réalisation des clichés, avec la méthode des éléments finis. Stockage dans snap2D/ ##
# Utilise fenics, éventuellement augmenté par dolfin #

#Caractérisation de la cellule élémentaire, tolérance ??

tol=1e-10

xinf=0.0
yinf=0.0
xsup=1.0
ysup=1.0
#c_x=0.5
#c_y=0.5

#determiner le domaine fixe pour interpoler la solution

class PeriodicBoundary(SubDomain):
 # Left boundary is "target domain" G
 def inside(self, x, on_boundary):
  return on_boundary and (near(x[0],xinf,tol) or near(x[1],yinf,tol))
 # Map right boundary (H) to left boundary (G)
 def map(self, x, y):
  if (near(x[0],xsup,tol)):
   y[0] = x[0] - 1.0
   y[1] = x[1]              
  else :
   y[0]=x[0]
   y[1] = x[1] - 1.0

res_fixe=80#résolution du maillage sans obstacle

domaine_fixe=Rectangle(Point(xinf,yinf),Point(xsup,ysup))
mesh_fixe=generate_mesh(domaine_fixe,res_fixe)
V_fixe=VectorFunctionSpace(mesh_fixe, "P", 2, constrained_domain=PeriodicBoundary())

#représentation graphique du maillage
plot(mesh_fixe)
plt.show()

#stockage dans Mesh2D

mesh_name='mesh_'+'fixe'+'_'+str(res_fixe)

rep_inc='Mesh2D'

########################################################################
### Commandes de Cyrille pour stocker et charger un maillage ###########
#ECRITURE                                                              #
#                                                                      #
#Name_mesh="%s/mesh_cylindre_fixe" %(repertoire_ecriture)+".xdmf"      #
#File(Name_mesh) << mesh                                               #
#                                                                      #
#LECTURE                                                               #
#                                                                      #
#Name_mesh="%s/mesh_cylindre_t=" %(repertoire_init)+str(T_init)+".xdmf"#
#mesh=Mesh(Name_mesh)                                                  #
########################################################################

### File(rep_inc+'/'+mesh_name+'.xdmf') << mesh_fixe #ne marche pas, contrairement à creation_fichier_pour_ecriture

def creation_fichier_pourecriture_champ_hdf5(repertoire_final,mesh):
	"Creation des dichiers pour ecrire la solution"
	ufile = File("%s/velocity.pvd" %(repertoire_final))                       
	USAVE=HDF5File(mesh.mpi_comm(),"%s/u_save.hdf5" %(repertoire_final), "w")
	return ufile,USAVE


ufile,USAVE=creation_fichier_pourecriture_champ_hdf5(rep_inc,mesh_fixe)
#file_rayon_ecriture = open("%s/rayon_ecriture.txt" %(repertoire_final), "w")
kfic=0


def ecriture_champ_hdf5(ufile,USAVE,u_n,kfic,file_rayon_ecriture,r):	
	"Ecriture des differents champs dans des fichiers"
	ufile << u_n
	USAVE.write(u_n,"khi",kfic)                                           
	file_rayon_ecriture.write(str(kfic)+"\t"+str(r)+"\n")
	return

#Créer un maillage : inclusion circulaire.




####################################

from DD_fun_obj import *



#tests pour les trois fonctions

c_x=0.01
c_y=0.01

r=0.35

res=25

mesh_c_r=creer_maill_circ([c_x,c_y],r,res)

plot(mesh_c_r)
plt.show()
plt.close()

figname='cx'+str(c_x)+'cy'+str(c_y)+'ray'+str(r)+'.png'
rep='Figures2D'
save_name=rep+'/'+figname

savefig(save_name)

#plt.show()

plt.close()


#définir : solveurEF
#définir








#Famille de cellules élémentaires : 8 clichés, inclusion circulaire

c_x=0.2
c_y=0.5

res=40

for i in range(1,3):#attention le rayon d'un cercle doit être non nul
 r=i*0.15
 mesh_c_r=creer_maill_circ([c_x,c_y],r,res)
 plot(mesh_c_r)
 # on enregistre ... comment ?
 plt.show()
 plt.close()
 # On pose et on résoud le problème aux éléments finis
 V=VectorFunctionSpace(mesh_c_r, 'P', 2, constrained_domain=PeriodicBoundary())
 ## On définit la bordure du domaine, sur laquelle intégrer le second membre "L" de l'équation en dimension finie
 l_cen=[]
 cen=[c_x,c_y]
 for i in range(-1,2):
  for j in range(-1,2):
   l_cen.append([cen[0]+i,cen[1]+j])
 ### Pour hériter de la méthode de marquage des facettes appartenant à Gamma_sf
 class inclusion_periodique(SubDomain):
  def inside(self,x,on_boundary):
   return any([between(sqrt((x[0]-c[0])**2+(x[1]-c[1])**2),(r-tol,r+tol)) for c in l_cen])
 ### Utilisation de la classe définie précédemment
 Gamma_sf = inclusion_periodique()
 boundaries = FacetFunction("size_t", mesh_c_r)
 boundaries.set_all(0)
 Gamma_sf.mark(boundaries, 5)
 ds = Measure("ds")(subdomain_data=boundaries)
 num_front_cercle=5 # numéro de la frontiere du cylindre
 ## On résoud le problème en fixant les condisions aux limites
 khi_bord=Constant((0., 0.))
 bc = DirichletBC(V, khi_bord, "x[0] < DOLFIN_EPS && x[1] < DOLFIN_EPS", "pointwise")
 normale = FacetNormal(mesh_c_r)
 nb_noeuds=V.dim()
 u = TrialFunction(V)
 v = TestFunction(V)
 a=tr(dot((grad(u)).T, grad(v)))*dx
 L=-dot(normale,v)*ds(num_front_cercle)
 ### Résolution
 u=Function(V)
 solve(a==L,u,bc)
 # Représentation graphique
 plot(u)
 plt.show()
 plt.close()

#Stockage



## Etape II : extrapolation des clichés, domaine_fixe ##

#Attention aux codes de Hassan : erreurs ... visibles sur les figures du rapport, s'il s'agit bien des snapshots extrapolés


## Etape III : en utilisant la méthode des snapshots, calcul de la POD et des coefficients aléatoires, toujours dans domaine_fixe ##


## Etape IV : Prédictions ##

# SE1 : ... #

## etc ##










