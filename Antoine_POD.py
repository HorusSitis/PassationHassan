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
	C += C.T
	
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
def mat_a_mat_phi(num_steps,U,C):

	egvl,egvct=np.linalg.eigh(C)
	valp=np.zeros((num_steps))
	A=np.zeros((num_steps,num_steps))
	for i in range(num_steps):
	    valp[i]=egvl[num_steps-1-i]
	    A[:,i]=egvct[:,num_steps-1-i]
	
	G=np.dot(U,A)
	return [valp,A,G]
