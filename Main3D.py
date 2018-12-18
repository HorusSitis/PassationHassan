# -*- coding: utf-8 -*-
### Une commande possible dans le terminal ###

#--- mpirun -np 8 python3 Main3D.py ---#
#--- affiche npfois 'pas encore' fait avec l'étape IV ---#

# Attention : on éxécute parallèlement 

### Code à lire : conditions ###

etape='EII'
# 'E0' à 'EIV' #
res_fixe=12
fixe_aff=False
res=12
Nsnap=8
rempUsnap='par8'#'seq'
c_x, c_y, c_z = 0.5, 0.5, 0.5
#r=0.35#pour une réalisation unique
npas_err=5

### Répertoire courant ###

#cd /home/amorea12/Documents/T_LaSIE/PassationHassan

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

from DDD_fun_obj import *
#import DDD_fun_obj as F3d
#from importlib import reload
#F3d=reload(F3d)

##############################################################################################################################
##############################################################################################################################
##############################################################################################################################

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

### Maillage sur le domaine Omega_fixe : aucune inclusion, porosité égale à 1 ###

res_fixe=6

domaine_fixe=Box(Point(xinf,yinf,zinf),Point(xsup,ysup,zsup))
mesh_fixe=generate_mesh(domaine_fixe,res_fixe)
V_fixe=VectorFunctionSpace(mesh_fixe, "P", 2, constrained_domain=PeriodicBoundary())

if fixe_aff==True :
 #représentation graphique du maillage
 plot(mesh_fixe)
 plt.show()
 plt.close()

if etape=='E0':
 #remplissage aléatoire'
 print('Remplissage fait')
### Etape I : réalisation des clichés, avec la méthode des éléments finis. Stockage dans snap2D/ ###
# Utilise fenics, éventuellement augmenté par dolfin #
elif etape=='EI':
 #
 #c_x=0.1
 #c_y=0.2
 #c_z=0.5
 #
 r=0.35
 #
 res=6
 #
 # Maillage raffiné unique #
 #
 mesh_c_r=creer_maill_sph([c_x,c_y,c_z],r,res)
 #
 plot(mesh_c_r)
 plt.show()
 plt.close()
 #
 #sys.exit()
 #
 # Champ khi unique #
 #
 #res=20
 #
 khi=snapshot_sph_per([c_x,c_y,c_z],r,res)
 #
 print(assemble(khi[1]*dx))
 print(assemble(khi[0]*dx))
 print(assemble(khi[2]*dx))
 print(khi((0.5,0.5,0.5)),khi((0.5,0.5,1.)),khi((0.5,0.5,0.0)))
 #
 plot(khi)
 plt.show()
 plt.close()
 #
 ## Erreur de périodicité ##
 #
 #
 print('Etape I- faite')
### Etape II : extrapolation des clichés, domaine_fixe ###
elif etape=='EII':
 #
 #c_x=0.5
 #c_y=0.5
 #c_z=0.5
 #
 #res=15
 #
 D_k=1.0
 #
 for i in range(1,2):
  r=i*0.35
  khi=snapshot_sph_per([c_x,c_y,c_z],r,res)
  ## chargement du snapshot pour l'indice courant
  # Extrapolation au domaine Omega_fixe : aucune inclusion, khi défini sur [0,1]times[0,1]
  khi.set_allow_extrapolation(True)
  khi_fixe=interpolate(khi,V_fixe)##rapide
  #u_fixe = project(u, V_fixe)##lent
  #print('erreur du gradient de khi :')
  #err_per_ind_01(grad(khi)[0,0],npas_err)##--> ne marche pas
  plot(grad(khi)[:,0])
  plt.show()
  plt.close()
  print('erreur sur la solution :')
  err_per_ind_01(khi,npas_err)
  plot(khi)
  #print('erreur sur la solution extrapolée :')
  #err_per_ind_01(khi_fixe,npas_err)
  #plot(khi_fixe)
  plt.show()
  plt.close()
 #
 ### Tenseur de diffusion homogénéisé
 ## Intégrale de khi sur le domaine fluide
 A=assemble(grad(khi)[0,0]*dx)
 B=assemble(grad(khi)[0,1]*dx)
 C=assemble(grad(khi)[0,2]*dx)
 E=assemble(grad(khi)[1,0]*dx)
 F=assemble(grad(khi)[1,1]*dx)
 G=assemble(grad(khi)[1,2]*dx)
 H=assemble(grad(khi)[2,0]*dx)
 K=assemble(grad(khi)[2,1]*dx)
 L=assemble(grad(khi)[2,2]*dx)
 Tkhi=np.array([[A,B,C],[E,F,G],[H,K,L]])
 ## Intégrale de l'identité sur le domaine fluide
 D=(1-4/3*pi*r**3)*np.eye(3)
 ## Calcul et affichage du tenseur Dhom
 Dhom=D_k*(D+Tkhi.T)
 #print(('Tenseur D_hom',Dhom))
 print(Dhom[0,0])
 print('Etape II- faite')
### Etape III : en utilisant la méthode des snapshots, calcul de la POD et des coefficients aléatoires, toujours dans domaine_fixe ###
elif etape=='EIII':
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
 print('Etape III- faite')
### Etape IV- : modèles réduits, comparaison avec la méthode des éléments finis ###
elif etape=='EIV':
 print('pas encore fait')


