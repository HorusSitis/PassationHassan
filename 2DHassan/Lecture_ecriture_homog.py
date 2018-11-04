# -*- coding: utf-8 -*-
"""
Created on Tue Jul 10 09:17:50 2018

@author: ghraieb
"""
import sys
from dolfin import *
import pylab as py

def ecriture_champ_hdf5(ufile,USAVE,u_n,kfic,file_rayon_ecriture,r):	
	"Ecriture des differents champs dans des fichiers"
	ufile << u_n
	USAVE.write(u_n,"Vitesse",kfic)                                           
	file_rayon_ecriture.write(str(kfic)+"\t"+str(r)+"\n")
	return


def creation_fichier_pourecriture_champ_hdf5(repertoire_final,mesh):
	"Creation des dichiers pour ecrire la solution"
	ufile = File("%s/velocity.pvd" %(repertoire_final))                       
	USAVE=HDF5File(mesh.mpi_comm(),"%s/u_save.hdf5" %(repertoire_final), "w")
	return ufile,USAVE


      
def lecture_fichiers_restart_hdf5(u_nmoins1,repertoire_init,T_init,mesh):      

	
	numfic,time = py.loadtxt("%s/rayon_ecriture.txt" %(repertoire_init), unpack=True)
	ilect=-1
	for i in range(len(time)):
#	    print(time[i],T_init)
	    if abs(time[i]-T_init)<1e-8 :
	      ilect=i
	if ilect==-1:
	    sys.exit("Erreur : Attention le temps de demarrage ne correspond pas a un temps de stockage \
			      -- voir dans le repertoire initial le fichier temps_ecriture.txt ")
      
	      
	
	URESTART=HDF5File(mesh.mpi_comm(),"%s/u_save.hdf5" %(repertoire_init), "r")


	dataset_urestart = "Vitesse/vector_%d"%ilect
	attr_urestart = URESTART.attributes(dataset_urestart)
	print 'Retrieving unmoins1 time step:', T_init, attr_urestart['timestamp'],dataset_urestart
	URESTART.read(u_nmoins1, dataset_urestart) 
	
	
		
	return



def lecture_solution_restart_hdf5(u_nmoins1,repertoire_init,T_init,mesh):      

	
	numfic,time = py.loadtxt("%s/rayon_ecriture.txt" %(repertoire_init), unpack=True)
	ilect=-1
	if abs(time-T_init)<1e-8 :
	      ilect=0
	if ilect==-1:
	    sys.exit("Erreur : Attention le temps de demarrage ne correspond pas a un temps de stockage \
			      -- voir dans le repertoire initial le fichier temps_ecriture.txt ")
      
	      
	
	URESTART=HDF5File(mesh.mpi_comm(),"%s/u_save.hdf5" %(repertoire_init), "r")


	dataset_urestart = "Vitesse/vector_%d"%ilect
	attr_urestart = URESTART.attributes(dataset_urestart)
	print 'Retrieving unmoins1 time step:', T_init, attr_urestart['timestamp'],dataset_urestart
	URESTART.read(u_nmoins1, dataset_urestart) 
	
	
		
	return