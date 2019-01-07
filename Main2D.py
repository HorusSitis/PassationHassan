






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
from importlib import reload

## Procédures pour la sauvegarde de fichiers ##

repertoire_parent="Res2D/"
#from LEc import *
import LEc as LE
#LE=reload(LE)

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

import DD_fun_obj as F2d
#F2d=reload(F2d)

### Codes éxécutés : cas d'une inclusion circulaire, le rayon du disque central est le paramètre pour la POD ###

#######################################################################################################################################
## Etape I : réalisation des clichés, avec la méthode des éléments finis. Calcul du tenseur d'homogénéisation. Stockage dans snap2D/ ##
#######################################################################################################################################

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

res_fixe=30#résolution du maillage sans obstacle

domaine_fixe=Rectangle(Point(xinf,yinf),Point(xsup,ysup))
mesh_fixe=generate_mesh(domaine_fixe,res_fixe)
V_fixe=VectorFunctionSpace(mesh_fixe, "P", 3, form_degree=1, constrained_domain=PeriodicBoundary())
###-------------------- Commandes pour l'écriture de fichiers, à déplacer dans le script éventuellement --------------------###

if [c_x,c_y]==[0.5,0.5]:
 suffixe="inc_centre/"
elif [c_x,c_y]==[0.0,0.0]:
 suffixe="coins/"

repertoire_final=repertoire_parent+suffixe

File(repertoire_parent+"mesh_circulaire.xml.gz") << mesh_fixe

kh_file,KH_SAVE=LE.creation_fichier_pourecriture_champ_hdf5(repertoire_final,mesh_fixe)

file_rayon_ecriture = open("%s/rayon_ecriture.txt" %(repertoire_final), "w")
kfic=1

###-------------------- pour des solutions interpolées --------------------###

#représentation graphique du maillage
plot(mesh_fixe)
plt.show()
plt.close()

# Famille de cellules élémentaires : 8 clichés, inclusion circulaire, paramétrée par le rayon du cercle

c_x=0.3
c_y=0.8

r=0.35

res=40

mesh_c_r=F2d.creer_maill_circ([c_x,c_y],r,res)

plot(mesh_c_r)
plt.show()
plt.close()

D_k=1.0

Npas=8

## Boucle pour la création des snapshots, avec un paramètre pouvant être le rayon d'une inclusion circulaire, ou l'emplacement de son centre ##
# Calcule aussi le tenseur de diffusion homogénéisé #

kfic=1
for i in range(7,1+Npas):#[0.111,0.211,0.316,0.423]:#,0.49]:#attention le rayon d'un cercle doit être non nul
 r=i*0.05
 c_x=0.3#5#3#i*0.1
 c_y=0.8#5#i*0.1
 khi_i=F2d.snapshot_circ_per([c_x,c_y],r,res)
 # Stockage des résultats avec un format hdf5
 ##LE.ecriture_champ_hdf5(kh_file,KH_SAVE,khi_i,kfic,file_rayon_ecriture,r,[c_x,c_y],res)
 print('Rayon :',r)
 #print(khi_i(0.8,0.0),khi_i(0.8,1.0))
 # Représentation graphique
 #plot(grad(khi_i)[:,0])
 #plt.show()
 #plt.close()
 for a in range(0,2):
  for b in range(0,2):
   p_gk=plot(-grad(khi_i)[a,b])
   plt.colorbar(p_gk)
   plt.show(p_gk)
   #plt.savefig("Figures2D/incccen_"+str(r)+"_dkhi"+str(a+1)+"_d"+str(b+1)+".png")
   plt.close()
 ##
 #plot(khi_i)
 #plt.show()
 #plt.close()
 for c in range(0,2):
  p_k=plot(khi_i[c])
  plt.colorbar(p_k)
  #plt.show(p_gk)
  plt.savefig("Figures2D/inc_c"+str(c_x)+str(c_y)+str(r)+"_khi"+str(c+1)+".png")
  plt.close()
 ##
 F2d.err_per_ind_01(khi_i,5)
 # Tenseur de diffusion homogénéisé
 ## Intégrale de khi sur le domaine fluide
 H=assemble(grad(khi_i)[0,0]*dx)
 M=assemble(grad(khi_i)[1,1]*dx)
 A=assemble(grad(khi_i)[0,1]*dx)
 C=assemble(grad(khi_i)[1,0]*dx)
 Tkhi=array([[H,A],[C,M]])
 ## Intégrale de l'identité sur le domaine fluide
 D=(1-pi*r**2)*np.eye(2)
 ## Calcul et affichage du tenseur Dhom
 Dhom=D_k*(D+Tkhi.T)
 #print(('Tenseur D_hom',Dhom))
 print(Dhom[0,0])
 # Stockage
 ## ...

#####################################################################################################################################
######################################### Etape II : extrapolation des clichés, domaine_fixe ########################################
#####################################################################################################################################

#----- Attention aux codes de Hassan : erreurs ... visibles sur les figures du rapport, s'il s'agit bien des snapshots extrapolés -----#
#----- Le même phénoène est observé avec mes codes, au 8-12-2018 -----#
#----- Peut-être s'agit-il d'effets dela géométrie de la cellule élémentaire : proxiité entre l'inclusion centrée et les bords de la cellule périodique -----#
#----- Le problème pourrait aussi venir de la méthode d'extrapolation -----#


#for i in range(1,1+Npas):
# r=i*0.2#05
# u=F2d.snapshot_circ_per([c_x,c_y],r,res)
# ## chargement du snapshot pour l'indice courant
# # Extrapolation au domaine Omega_fixe : aucune inclusion, khi défini sur [0,1]times[0,1]
# u.set_allow_extrapolation(True)
# u_fixe=interpolate(u,V_fixe)##rapide
# #u_fixe = project(u, V_fixe)##lent
# plot(u_fixe)
# plt.show()
# plt.close()

## On recharge les clichés stockés à l'étape I, aec le format hdf5
if [c_x,c_y]==[0.5,0.5]:
 suffixe="inc_centre/"
elif [c_x,c_y]==[0.0,0.0]:
 suffixe="coins/"

repertoire_lecture=repertoire_parent+suffixe

Npas=8
#c_x,c_y=0.5,0.5
res=40
rayons=np.arange(0.05,0.45,0.05)
rayons=np.linspace(0.05,0.40,8)

for i in range(1,1+Npas):
 r=rayons[i]
 # Création de l'espace dans lequel vit le cliché à charger
 mesh_c_r=creer_maill_circ([c_x,c_y],r,res)
 V=VectorFunctionSpace(mesh_c_r, 'P', 2, constrained_domain=PeriodicBoundary())
 khi_i=Function(V)
 # Chargement du cliché, réalisé à l'étape I
 LE.lecture_fichiers_restart_hdf5(khi_i,repertoire_lecture,r,mesh_c_r)


####################################################################################################################################
## Etape III : en utilisant la méthode des snapshots, calcul de la POD et des coefficients aléatoires, toujours dans domaine_fixe ##
####################################################################################################################################

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

vp_A_phi=pod.mat_a_mat_phi(Nsnap,Usnap,C,V_fixe,'n2')
vp_A_phi=pod.mat_a_mat_phi(Nsnap,Usnap,C,V_fixe,'L2')
#vp_A_phi=pod.mat_a_mat_phi(Nsnap,Usnap,C,'')

val_propres=vp_A_phi[0]
Aleat=vp_A_phi[1]
## Attention les objets rangés dans tableau suivant sont des vecteurs
Phi_prime_v=vp_A_phi[2]

## Tests : orthogonalité ou orthonrmalité de Phi_prime
ui=Function(V_fixe)
uj=Function(V_fixe)

## Orthogonalité
for i in range(Nsnap-1):
 ui.vector().set_local(Phi_prime_v[:,i])
 for j in range(i+1,Nsnap):
  uj.vector().set_local(Phi_prime_v[:,j])
  scal=assemble(dot(ui,uj)*dx)
  print(scal)

## Norme des vacteurs dela base POD, L2 ou n2
for i in range(Nsnap):
 ui.vector().set_local(Phi_prime_v[:,i])
 scal=assemble(dot(ui,ui)*dx)
 norme_L2=sqrt(scal)
 ###
 norme_q=0
 l=len(Phi_prime_v[:,i])
 for k in range(l):
  norme_q+=Phi_prime_v[k,i]**2
 norme_2=sqrt(norme_q)
 #print('norme 2 :',norme_q)
 print('norme L2 :',norme_L2)
 #print('quotient n2/L2 :',scal/norme_q)

# Représentation graphique des vecteurs de POD :

## Type de données : on veut calculer les fonctions phi_prime_i 
## Représentation graphique des phi_prime_i :

phi=Function(V_fixe)
for i in range(Nsnap):
 phi.vector().set_local(Phi_prime_v[:,i])
 plot(phi)
 plt.show()
 plt.close()

# Energie et énergie cumulée des modes spatiaux, choix du nombre de modes

## Energie et énergie cumulée, avec les valeurs propres de la matrice de corrélation temporelle
ener_pour=pod.energie_pourcentage(val_propres)[0]
ener_pour_cumul=pod.energie_pourcentage(val_propres)[1]

absc=np.arange(1,Nsnap+1,1)

plt.plot(absc,ener_pour)
plt.xlabel('valeurs propres')
plt.ylabel('pourcentage_energie')
plt.show()

plt.plot(absc,ener_pour_cumul)
plt.xlabel('valeurs propres')
plt.ylabel('pourcentage_energie_cumule')
plt.show()

## Choix du nombre de modes, avec une valeur seuil d'énergie à atteindre avec les vacteurs de la base POD
nb_modes=0

seuil_ener=99.999

i=0
while ener_pour_cumul[i]<seuil_ener:
 nb_modes=i+1
 i+=1

### 8 snapshots : 4 modes pour un seuil de 99.9%
### 8 snapshots : 7 modes pour un seuil de 99,999%

#################################################################################################
## Etape IV : Prédictions. Choisir les paramètres du problème à résoudre par le modèle réduit. ##
#################################################################################################

c_x=
c_y=
r=

# SE1 : projection de la base POD sur le nouveau domaine #





## etc ##










