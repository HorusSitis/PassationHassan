# -*- coding: utf-8 -*-
"""
Created on Thu Jun  7 12:20:45 2018

@author: ghraieb
"""

from fenics import *
from math import *
import numpy as np
import matplotlib.pyplot as plt

tol=1e-10

xinf=0.0
yinf=0.0
zinf=0.0
xsup=1.0
ysup=1.0
zsup=1.0

#########################################
####Matrice de Corrélation temporelle####
#########################################

def mat_corr_temp(V,num_steps,U):
	
	C=np.zeros((num_steps,num_steps),dtype=np.float64)
	S=np.zeros((num_steps))
	ui=Function(V)
	uj=Function(V)

	for k in range(num_steps-1):
		ui.vector().set_local(U[:,k])
		for i in range(k+1,num_steps):
			uj.vector().set_local(U[:,i])
			C[k,i]=assemble(dot(ui,uj)*dx)         
	C=C*(1.0/num_steps)
	C += C.T## on a effectué seulement num_steps*(num_step+1)/2 produits entre les fonctions ui-j
	
	for k in range(num_steps):
		ui.vector().set_local(U[:,k])
		S[k]=assemble(dot(ui,ui)*dx)
	S=S*(1.0/num_steps)
	N=np.diag(S)

	return C+N

##########################################################################################
#######################Matrice des coefficients aléatoires A (ai(t))######################
####Matrice de la base totale sans prendre en compte l'énergie de chaque vecteur (phi)####
##########################################################################################

def mat_a_mat_phi(num_steps,U,C,V,base_POD_normee):
 ## Résolution du problème aux valeurs-vecteurs propres, avec la matrice des corrélations temporelles
 egvl,egvct=np.linalg.eigh(C)
 valp=np.zeros((num_steps))
 A=np.zeros((num_steps,num_steps))
 ## Remplissage pour les sorties A et vp
 for i in range(num_steps):
  valp[i]=egvl[num_steps-1-i]
  A[:,i]=egvct[:,num_steps-1-i]
 ## Calcul des vaceurs POD avec la matrice U des snapshots et l'estimation A
 Phi_prime=np.dot(U,A)
 ### On norme la base POD, optionnel
 phi=Function(V)
 if base_POD_normee=='n_2' :
  for i in range(num_steps):
   phi_prime_i=Phi_prime[:,i]
   norme_q=0
   l=len(phi_prime_i)
   for k in range(l):
    norme_q=norme_q+phi_prime_i[k]**2
   norme_2=sqrt(norme_q)
   print(norme_2)
   phi_prime_i=phi_prime_i/norme_2
 elif base_POD_normee=='L2':
  for i in range(num_steps):
   phi.vector().set_local(Phi_prime[:,i])
   scal=assemble(dot(phi,phi)*dx)
   norme_L2=sqrt(scal)
   Phi_prime[:,i]=Phi_prime[:,i]*(1.0/norme_L2)
 ## Résultats : valeurs propres ... ordre ? ; matrice de coefficients aléatoires A et base POD Phi_prime
 return [valp,A,Phi_prime]

##########################################################################################
##################### Energie et énergie cumulée des valeurs propres #####################
##########################################################################################

def energie_pourcentage(vp):
 R_dim=len(vp)
 s_t=0
 ener_pour=np.zeros((R_dim))
 ener_pour_cumul=np.zeros((R_dim))
 for k in range(R_dim):
  s_t=s_t+vp[k]
 for i in range(R_dim):
  s=vp[i]
  s_n=0
  for j in range(i+1):
   s_n=s_n+vp[j]
  ener_pour[i]=(s/s_t)*100
  ener_pour_cumul[i]=(s_n/s_t)*100
 return([ener_pour,ener_pour_cumul])

##################################################################################################################################################################
##################### Calcul des coefficients du modèle réduit, étant donné le domaine de définition du nouveau champ de vecteurs à calculer #####################
##################################################################################################################################################################

#V_nouv
#Phi_prime_v


tol=1e-10
def calc_Ab_2D(V_nouv,mesh_nouv,Phi_nouv_v,r_nouv,cen,nb_modes):
 A=np.zeros((nb_modes,nb_modes))
 b=np.zeros(nb_modes)
 ### Fonctions à définir pour calculer les coefficients des deux tenseurs, qui dépendent de la métrique de l'espace des fonctions test
 phi_nouv_k=Function(V_nouv)
 phi_nouv_i=Function(V_nouv)
 # boucle pour le calcul de la matrice de coefficients
 for k in range(nb_modes):
  phi_nouv_k.vector().set_local(Phi_nouv_v[:,k])
  for i in range(nb_modes):
   phi_nouv_i.vector().set_local(Phi_nouv_v[:,i])
   # On calcule le coefficient Aki
   A[k,i]=assemble(tr(dot((grad(phi_nouv_k)).T, grad(phi_nouv_i)))*dx)
 # création de l'interface solide-fluide
 l_cen=[]
 for i in range(-1,2):
  for j in range(-1,2):
   l_cen.append([cen[0]+i,cen[1]+j])
 r=r_nouv
 class inclusion_periodique(SubDomain):
  def inside(self,x,on_boundary):
   return (on_boundary and any([between((x[0]-c[0]), (-r-tol, r+tol)) for c in l_cen]) and any([between((x[1]-c[1]), (-r-tol, r+tol)) for c in l_cen]))
 Gamma_sf=inclusion_periodique()
 boundaries = MeshFunction("size_t", mesh_nouv, mesh_nouv.topology().dim()-1)
 boundaries.set_all(1)
 Gamma_sf.mark(boundaries, 7)
 ds = Measure("ds")(subdomain_data=boundaries)
 num_ff=1
 num_front_inc=7
 normale=FacetNormal(mesh_nouv)
 # boucle pour le calcul du second membre du problème linéaire MOR
 for i in range(nb_modes):
  phi_nouv_i.vector().set_local(Phi_nouv_v[:,i])
  b[i]=assemble(dot(normale,phi_nouv_i)*ds(num_front_inc))
 return([A,b])

def calc_Ab_3D(V_nouv,mesh_nouv,Phi_nouv_v,r_nouv,origin,nb_modes,config):
 A=np.zeros((nb_modes,nb_modes))
 b=np.zeros(nb_modes)
 ### Fonctions à définir pour calculer les coefficients des deux tenseurs, qui dépendent de la métrique de l'espace des fonctions test
 phi_nouv_k=Function(V_nouv)
 phi_nouv_i=Function(V_nouv)
 # boucle pour le calcul de la matrice de coefficients
 for k in range(nb_modes):
  phi_nouv_k.vector().set_local(Phi_nouv_v[:,k])
  for i in range(nb_modes):
   phi_nouv_i.vector().set_local(Phi_nouv_v[:,i])
   # On calcule le coefficient Aki
   A[k,i]=assemble(tr(dot((grad(phi_nouv_k)).T, grad(phi_nouv_i)))*dx)
 # création de l'interface solide-fluide
 r=r_nouv
 if config=='sph_un':
  l_cen=[origin]
  print(l_cen)
  class inclusion_periodique(SubDomain):
   def inside(self,x,on_boundary):
    return (on_boundary and any([between((x[0]-c[0]), (-r-tol, r+tol)) for c in l_cen]) and any([between((x[1]-c[1]), (-r-tol, r+tol)) for c in l_cen]) and any([between((x[2]-c[2]), (-r-tol, r+tol)) for c in l_cen]))
 elif config=='cyl_un':
  l_axe=[origin]
  print(l_axe)
  class inclusion_periodique(SubDomain):
   def inside(self,x,on_boundary):
    return (on_boundary and any([between((x[0]-c[0]), (-r-tol, r+tol)) for c in l_axe]) and any([between((x[2]-c[2]), (-r-tol, r+tol)) for c in l_axe]))
 Gamma_sf=inclusion_periodique()
 boundaries = MeshFunction("size_t", mesh_nouv, mesh_nouv.topology().dim()-1)
 boundaries.set_all(1)
 Gamma_sf.mark(boundaries, 7)
 ds = Measure("ds")(subdomain_data=boundaries)
 num_ff=1
 num_front_inc=7
 normale=FacetNormal(mesh_nouv)
 # boucle pour le calcul du second membre du problème linéaire MOR
 for i in range(nb_modes):
  phi_nouv_i.vector().set_local(Phi_nouv_v[:,i])
  b[i]=assemble(dot(normale,phi_nouv_i)*ds(num_front_inc))
 return([A,b])

### Valable en toutes dimensions ###

def calc_Ab_compl(V_nouv,mesh_nouv,Phi_nouv_v,nb_modes,test_snap):
 A=np.zeros((nb_modes,nb_modes))
 b=np.zeros(nb_modes)
 ## Fonctions à définir pour calculer les coefficients des deux tenseurs, qui dépendent de la métrique de l'espace des fonctions test
 phi_nouv_k=Function(V_nouv)
 phi_nouv_i=Function(V_nouv)
 ## Boucle pour le calcul de la matrice de coefficients
 for k in range(nb_modes):
  phi_nouv_k.vector().set_local(Phi_nouv_v[:,k])
  for i in range(nb_modes):
   phi_nouv_i.vector().set_local(Phi_nouv_v[:,i])
   # On calcule le coefficient Aki
   A[k,i]=assemble(tr(dot((grad(phi_nouv_k)).T, grad(phi_nouv_i)))*dx)
 ## On définit la bordure du domaine, sur laquelle intégrer le second membre "L" de l'équation en dimension finie
 boundaries = MeshFunction("size_t", mesh_nouv, mesh_nouv.topology().dim()-1)
 #boundaries = MeshFunction('size_t', mesh, mesh_name+"_facet_region"+".xml")
 ds = Measure("ds")(subdomain_data=boundaries)
 ## Marquage des bordures pour la condition de Neumann
 if test_snap=='solid_1':
  num_front_inc=1
  class SolidBoundary(SubDomain):
   def inside(self, x, on_boundary):
    return on_boundary and not(near(x[0],xinf,tol) or near(x[0],xsup,tol) or near(x[1],yinf,tol) or near(x[1],ysup,tol))
  Gamma_sf = SolidBoundary()
  print('Gamma sf ne coupe pas le bord du carré')
  boundaries.set_all(0)
  Gamma_sf.mark(boundaries, 1)
 elif test_snap=='solid_2':
  #
  num_front_inc=11
  print('Gamma sf coupe le bord du carré')
  #
 ## On intègre les vecteurs POD pour obtenir les coefficients du modèle réduit
 normale=FacetNormal(mesh_nouv)
 # boucle pour le calcul du second membre du problème linéaire MOR
 for i in range(nb_modes):
  phi_nouv_i.vector().set_local(Phi_nouv_v[:,i])
  b[i]=assemble(dot(normale,phi_nouv_i)*ds(num_front_inc))
 return([A,b])

## Dimension 3 : même algorithme ##

def calc_Ab_compl_3D(V_nouv,mesh_nouv,Phi_nouv_v,nb_modes):#,test_snap):
 A=np.zeros((nb_modes,nb_modes))
 b=np.zeros(nb_modes)
 ## Fonctions à définir pour calculer les coefficients des deux tenseurs, qui dépendent de la métrique de l'espace des fonctions test
 phi_nouv_k=Function(V_nouv)
 phi_nouv_i=Function(V_nouv)
 ## Boucle pour le calcul de la matrice de coefficients
 for k in range(nb_modes):
  phi_nouv_k.vector().set_local(Phi_nouv_v[:,k])
  for i in range(nb_modes):
   phi_nouv_i.vector().set_local(Phi_nouv_v[:,i])
   # On calcule le coefficient Aki
   A[k,i]=assemble(tr(dot((grad(phi_nouv_k)).T, grad(phi_nouv_i)))*dx)
 ## On définit la bordure du domaine, sur laquelle intégrer le second membre "L" de l'équation en dimension finie
 boundaries = MeshFunction("size_t", mesh_nouv, mesh_nouv.topology().dim()-1)
 #boundaries = MeshFunction('size_t', mesh, mesh_name+"_facet_region"+".xml")
 ds = Measure("ds")(subdomain_data=boundaries)
 ## Marquage des bordures pour la condition de Neumann
 num_front_inc=1
 class SolidBoundary(SubDomain):
  def inside(self, x, on_boundary):
   return on_boundary and not(near(x[0],xinf,tol) or near(x[0],xsup,tol) or near(x[1],yinf,tol) or near(x[1],ysup,tol) or near(x[2],zinf,tol) or near(x[2],zsup,tol))
 Gamma_sf = SolidBoundary()
 #print('Gamma sf ne coupe pas le bord du cube')
 boundaries.set_all(0)
 Gamma_sf.mark(boundaries, 1)
 ## On intègre les vecteurs POD pour obtenir les coefficients du modèle réduit
 normale=FacetNormal(mesh_nouv)
 # boucle pour le calcul du second membre du problème linéaire MOR
 for i in range(nb_modes):
  phi_nouv_i.vector().set_local(Phi_nouv_v[:,i])
  b[i]=assemble(dot(normale,phi_nouv_i)*ds(num_front_inc))
 return([A,b])

