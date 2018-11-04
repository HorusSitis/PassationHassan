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
	return valp,A,G



#######################################################################
##########Choix du nb de mode à (10^4) près de lamda maximale##########
#######################################################################
def calc_nb_modes(num_steps,eigenvalues):

	maximum=eigenvalues[0];

	n=1;
	for i in range(num_steps-1):
		if ((eigenvalues[i+1])>(maximum*pow(10,(-4)))):
			n=n+1
		else :
			break
	
	return n



#########################################################
####Test de l'orthogonalité de la famille de base phi####
#########################################################
def test_orthog(V,G,num_steps):

	ui=Function(V)
	uj=Function(V)
   ##################################################
   ##########Pour orthonormer la base################	
   ##################################################
	
#	for i in range(num_steps):
#	       ui.vector().set_local(G[:,i])
#	       scal=assemble(dot(ui,ui)*dx)
#	       G[:,i]=G[:,i]/sqrt(scal)

	for i in range(num_steps):
	        ui.vector().set_local(G[:,i])
	        for j in range(num_steps):
	                uj.vector().set_local(G[:,j])
	                scal=assemble(dot(ui,uj)*dx)
	                if (i<2):
                         print(i,j,scal)

	
	return



#########################################################
#########calcul de l'erreur en champ instantanné#########
#########################################################
def erreur_nbmode_inst(V,U,nbmode,num_steps,nb_noeuds,A,G):

	erreur=np.zeros((nbmode,num_steps))
	Urec=np.zeros((nb_noeuds,num_steps))
	
	ui=Function(V)
	vi=Function(V)	
	for j in range(nbmode):
	    
	    for i in range(num_steps):
	        a=A[i,0:j+1]
	        phi=G[:,0:j+1]
	        Urec[:,i]=np.dot(phi,a)
	        numerateur=U[:,i]-Urec[:,i]
	        denominateur=U[:,i]
	        ui.vector().set_local(numerateur)
	        vi.vector().set_local(denominateur)
	        modul_ui_carree=assemble(dot(ui,ui)*dx)
	        modul_vi_carree=assemble(dot(vi,vi)*dx)
	        erreur[j,i]=np.sqrt(modul_ui_carree)/np.sqrt(modul_vi_carree)
	
	errnbm=np.zeros((nbmode))
	for i in range(nbmode):
	    som=0.0
	    for j in range(num_steps):
	        som=som+erreur[i,j]
	    errnbm[i]=som/(num_steps)

####pour afficher les valeurs des erreurs####
#	    print("erreur du mode ",i,errnbm[i])              
	
	
	return errnbm



#########################################################
##########calcul de l'erreur en champ fluctuant##########
#########################################################
def erreur_nbmode(V,Ubar,Utab,nbmode,num_steps,nb_noeuds,A,G):

	erreur=np.zeros((nbmode,num_steps))
	Uprim_rec=np.zeros((nb_noeuds,num_steps))
	Urec=np.zeros((nb_noeuds,num_steps))
	
	ui=Function(V)
	vi=Function(V)	
	for j in range(nbmode):
	    
	    for i in range(num_steps):
	        a=A[i,0:j+1]
	        phi=G[:,0:j+1]
	        Uprim_rec[:,i]=np.dot(phi,a)
	        Urec[:,i]=Uprim_rec[:,i]+Ubar
	        numer=Utab[:,i]-Urec[:,i]
	        denom=Utab[:,i]
	        ui.vector().set_local(numer)
	        vi.vector().set_local(denom)
	        modul_ui_carree=assemble(dot(ui,ui)*dx)
	        modul_vi_carree=assemble(dot(vi,vi)*dx)
	        erreur[j,i]=np.sqrt(sous)/np.sqrt(sous1)
	
	errnbm=np.zeros((nbmode))
	for i in range(nbmode):
	    som=0
	    for j in range(num_steps):
	        som=som+erreur[i,j]
	    errnbm[i]=som/(num_steps)

####pour afficher les valeurs des erreurs####
#	    print("erreur du mode ",i,errnbm[i])              
	
	
	return errnbm