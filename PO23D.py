# -*- coding: utf-8 -*-
"""
Created on Thu Jun  7 12:20:45 2018

@author: ghraieb
"""

from fenics import *
from math import *
import numpy as np
import matplotlib.pyplot as plt

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






