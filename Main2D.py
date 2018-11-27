






###############################################################################################################
############################### Génération de structures périodiques aléatoires ###############################
###############################################################################################################

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

##500 ; tps 30 000 ; 4600 inclusions : environ 1h sur MECALAC 8 x 2,90 GHz
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

### Répertoire courant ###

cd /home/amorea12/Documents/T_LaSIE/PassationHassan

### Paquets à importer ###

import numpy as np
#import random as rd
import pylab as pl
from pylab import *
#import os
import multiprocessing

## Paquets spécifiques à POD-MOR ##

from fenics import *
from dolfin import *
from mshr import *
import matplotlib.pyplot as plt
from math import sqrt
import sys

#Imports divers : fonctions, objets etc

#from DD_fun_obj import *#attention pas de chiffre au début du nom d'un paquet
#r=0.0
## problème à l'import : r prend la valeur index(6)

from importlib import reload

import DD_fun_obj as F2d

F2d=reload(F2d)

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
   # pour le coin supérieur droit : cohérence avec x[1] == 1.0
   if (near(x[1],ysup,tol)):
    y[0] = x[0] - 1.0
    y[1] = x[1] - 1.0
   # pour le reste du bord droit
   else:
    y[0] = x[0] - 1.0
    y[1] = x[1]   
  elif            
  else :
   y[0]=x[0]
   y[1] = x[1] - 1.0

res_fixe=30#résolution du maillage sans obstacle

domaine_fixe=Rectangle(Point(xinf,yinf),Point(xsup,ysup))
mesh_fixe=generate_mesh(domaine_fixe,res_fixe)
V_fixe=VectorFunctionSpace(mesh_fixe, "P", 2, constrained_domain=PeriodicBoundary())

#représentation graphique du maillage
plot(mesh_fixe)
plt.show()
plt.close()

#stockage dans Mesh2D

mesh_name='mesh_'+'fixe'+'_'+str(res_fixe)

rep_inc='Mesh2D'
#sys.exit("fin temporaire")

# Famille de cellules élémentaires : 8 clichés, inclusion circulaire, paramétrée par le rayon du cercle

c_x=0.5
c_y=0.5

res=25

D_k=1.0

Npas=4

#l_err_abs=[]
#l_err_rel=[]

for i in range(1,6):#[0.111,0.211,0.316,0.423]:#,0.49]:#range(1,2):#9):#attention le rayon d'un cercle doit être non nul
 r=i*0.05
 c_x=0.5#i*0.1
 c_y=0.5#i*0.1
 u=F2d.snapshot_circ_per([c_x,c_y],r,res)
 # Représentation graphique
 #plot(u)
 #plt.show()
 #plt.close()
 # erreurs de périodicité, point par point
 pas=1/Npas
 #l_x_0=array([u((0,k*pas)) for k in range(0,Npas)])
 #l_x_1=array([u((1,k*pas)) for k in range(0,Npas)])
 #list_per_x=[sum((u((1,pas*k))-u((0,pas*k)))**2) for k in range(0,Npas)]
 #l_x_diff=l_x_1-l_x_0
 #El2_per_x=sqrt(sum(list_per_x)/Npas)/1.0#den[0]
 #list_0_x=[sum(u((0,pas*k))**2) for k in range(0,Npas)]
 #El2_rel_per_x=sqrt(sum(list_per_x)/Npas)/sqrt(sum(list_0_x)/Npas)
 #print('Différence en x :',l_x_diff)
 #print('Différence et valeur quadratiques en x :',list_per_x,list_0_x)
 #print('Différence en moyenne quadratique pour x :',El2_per_x)
 #print('Différence relative en moyenne quadratique pour x :',El2_rel_per_x)
 #l_y_0=array([u((0.1*k,0)) for k in range(0,10)])
 #l_y_1=array([u((0.1*k,1)) for k in range(0,10)])
 #l_y_diff=l_y_1-l_y_0
 #print('Différence en y :',l_y_diff)
 ##
 ### Extrapolation au domaine entier : [0,1]^2 ###
 #u_fixe=u
 ## erreurs d'extrapolation, absolue et relative
 u.set_allow_extrapolation(True)
 u_fixe=interpolate(u,V_fixe)
 #plot(u_fixe)
 #plt.show()
 #plt.close()
 # Périodicité
 print('khiprime au bord vertical :')
 for k in range(0,1+Npas):
  print(u_fixe((0,pas*k)),u_fixe((1,pas*k)))
 print('khiprime au bord horizontal :')
 for k in range(0,1+Npas):
  print(u_fixe((pas*k,0)),u_fixe((pas*k,1)))
 #print('u',err_per_01(u,'l2',100,'rel'))
 #print('ufixe',err_per_01(u_fixe,'l2',100,'rel'))
 # Erreur d'interpolation
 u_mat=np.zeros((1+Npas,1+Npas))
 u_diff=np.zeros((1+Npas,1+Npas))
 #for k in range(0,1+Npas):
 # for l in range(0,1+Npas):
 #  u_diff[k,l]=sum([x**2 for x in u((pas*k,pas*l))])-sum([x**2 for x in u_fixe((pas*k,pas*l))])
 #  u_mat[k,l]=sum([x**2 for x in u((pas*k,pas*l))])
 #u_arr=np.reshape(u_mat,(1+Npas)*(1+Npas))
 #u_diff_ligne=np.reshape(u_diff,(1+Npas)*(1+Npas))
 #u_rel_ligne=u_diff_ligne/u_arr
 #diff_unif=max(u_diff_ligne)
 #err_u_rel=max(u_rel_ligne)#diff_unif/u_min_unif
 #l_err_abs.append(('rayon',r,'abs',diff_unif))
 #l_err_rel.append(('rayon',r,'erreur relative maximum',err_u_rel))
 #print('erreur interpolation',diff_unif)
 #print('erreur relative interpolation',err_u_rel)
 #print('centre de la cellule',u((0.5,0.5)),u_fixe((0.5,0.5)))
 ##
 # Tenseur d'homogéisation
 ## Intégrale de khi sur le domaine fluide
 H=assemble(grad(u)[0,0]*dx)
 M=assemble(grad(u)[1,1]*dx)
 A=assemble(grad(u)[0,1]*dx)
 C=assemble(grad(u)[1,0]*dx)
 Tkhi=array([[H,A],[C,M]])
 ## Intégrale de l'identité sur le domaine fluide
 D=(1-pi*r**2)*np.eye(2)
 ## Calcul et affichage du tenseur Dhom
 Dhom=D_k*(D+Tkhi.T)
 #print(('Tenseur D_hom',Dhom))
 print(Dhom[0,0])
 # Stockage
 ## ...

## Erreurs d'interpolation : c_x, c_y = 0.5,0.5
#print('abs',l_err_abs)
for i in range(0,8):
 print(l_err_rel[i])
## 0.05 à 0.4, suite arithmétique :
#('rayon', 0.05, 'erreur relative maximum', 0.0011250471700128834)
#('rayon', 0.1, 'erreur relative maximum', 0.0024562097480711302)
#('rayon', 0.15000000000000002, 'erreur relative maximum', 0.0079380966355682512)
#('rayon', 0.2, 'erreur relative maximum', 0.031271739691666682)
#('rayon', 0.25, 'erreur relative maximum', 0.021105098011912495)
#('rayon', 0.30000000000000004, 'erreur relative maximum', 0.032643098011902655)
#('rayon', 0.35000000000000003, 'erreur relative maximum', 0.059268219555300893)
#('rayon', 0.4, 'erreur relative maximum', 0.16544547860510037)
## -- End pasted text --

## Erreurs d'interpolation : c_x, c_y = 0.2,0.2 ; suite arithmétique de raison etpremier terme 0.05
#print('abs',l_err_abs)
for i in range(0,5):
 print(l_err_rel[i])
#('rayon', 0.05, 'erreur relative maximum', 0.30378848343770709)
#('rayon', 0.1, 'erreur relative maximum', 0.49584821840470178)
#('rayon', 0.15000000000000002, 'erreur relative maximum', 0.57479441504707507)
#('rayon', 0.2, 'erreur relative maximum', 1.4487937386039291e-09)
#('rayon', 0.25, 'erreur relative maximum', 3.0524331809039409e-08)
## -- End pasted text --

## Etape II : extrapolation des clichés, domaine_fixe ##

# Attention aux codes de Hassan : erreurs ... visibles sur les figures du rapport, s'il s'agit bien des snapshots extrapolés

c_x=0.5
c_y=0.5

res=25

for i in range(1,2):
 r=i*0.2#05
 u=F2d.snapshot_circ_per([c_x,c_y],r,res)
 ## chargement du snapshot pour l'indice courant
 # Extrapolation au domaine Omega_fixe : aucune inclusion, khi défini sur [0,1]times[0,1]
 u.set_allow_extrapolation(True)
 u_fixe=interpolate(u,V_fixe)##rapide
 ## Erreur de périodicité du snapshot calculé par éléments finis
 print('Erreur l2 absolue :',F2d.err_per_01(u,'l2',100,''))
 print('Erreur l2 relative :',F2d.err_per_01(u,'l2',100,'rel'))
 print('Erreur infty absolue :',F2d.err_per_01(u,'infty',100,''))
 print('Erreur infty relative :',F2d.err_per_01(u,'infty',100,'rel'))
 #u_fixe = project(u, V_fixe)##lent
 plot(u_fixe)
 plt.show()
 plt.close()

## Etape III : en utilisant la méthode des snapshots, calcul de la POD et des coefficients aléatoires, toujours dans domaine_fixe ##

#mesh_fixe = Mesh("Solutions_homog_interp_circulaire/mesh_circulaire.xml.gz")

import PO23D as pod
pod = reload(pod)
#from Antoine_POD import *

# Domaine d'interpolation et matrice des snapshots

c_x=0.5
c_y=0.5

res=25

Nsnap=8

#V_fixe=VectorFunctionSpace(mesh_fixe, "P", 2, constrained_domain=PeriodicBoundary())
nb_noeuds = V_fixe.dim()

Usnap=np.zeros((nb_noeuds,Nsnap))

def inter_snap_ray(n):
 r=n*0.05
 # Cliché sur le domaine avec inclusion
 u=F2d.snapshot_circ_per([c_x,c_y],r,res)
 # Extrapolation du cliché : khi prime
 u.set_allow_extrapolation(True)
 u_fixe=interpolate(u,V_fixe)
 # Forme vectorielle de la solution EF
 u_fixe_v=u_fixe.vector().get_local()
 return([n,u_fixe_v])#[n,u_fixe])

# Génération séquentielle des snapshots

for n in range(1,1+Nsnap):
 u_fixe_v=inter_snap_ray(n)[1]
 # Remplissage de la matrice des snapshots
 Usnap[:,n-1]=u_fixe_v#.vector().get_local()

## UsnapSeq=Usnap

# Génération parallèle des snapshots

pool=multiprocessing.Pool(processes=8)

Uf_par=pool.map(inter_snap_ray,(n for n in range(1,1+Nsnap)))

for n in range(1,1+Nsnap):
 for i in range(0,Nsnap):
  if Uf_par[i][0]==n:
   u_fixe_v=Uf_par[i][1]
   Usnap[:,n-1]=u_fixe_v#.vector().get_local()

## UsnapPar=Usnap

# matrice de corrélation

## Usnap=UsnapSeq
## Usnap=UsnapPar

C=pod.mat_corr_temp(V_fixe,Nsnap,Usnap)

# Calcul des coefficients aléatoires et la base POD

vp_A_phi=pod.mat_a_mat_phi(Nsnap,Usnap,C)

val_propres=vp_A_phi[0]
Aleat=vp_A_phi[1]
phi=vp_A_phi[2]



## Etape IV : Prédictions ##

# SE1 : ... #

## etc ##










