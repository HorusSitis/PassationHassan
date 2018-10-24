### Encodage ###


### Paquets à importer ###

import numpy as np
import random as rd

#import pylab as pl
#from pylab import *


#import os


## Remplissage avec des sphères
#Entrées : une norme, une distance minimu entre les inclusions, les lois de probabilité pour les rayons des inclusions et pour chaque phase, les paramètres choisis pour ces lois, une fraction volumique pour chaque phase incluse dans l'espace ambiant initialement vide et un triplet pour le domaine spatial étudié
#Sortie : liste unique, des rayons et centres des inclusions, indexés par le numéro de phase.

#Principe : A chaque tour de boucle, on crée une sphère pour une phase donnée, et on vérifie que son centre est suffisamment loin de tous les autres centres déjà créés.
#Cas périodique ?

##Remarque : l'arrêt de la procédure principale n'est pas prouvé

#A=np.zeros(27)
#A=np.reshape(A,(3,3,3))

#Etape préliminaire : fonctions pour générer des nombres aléatoires dans des boucles. Lois éventuellement tronquées pour obtenir des valeurs positives, on cherche à simuler des rayons

def g_norm(par):return(rd.gauss(par[0],par[1]))

#Géométrie choisie avec la première variable : formules pour la distance et le volume

##Exemple : geom=[eucl3D,volB] liste de deux fonctions à appeler dans l'algorithme

def eucl2D(a,b):return(np.sqrt((a[0]-b[0])**2+(a[1]-b[1])**2))

def vol2D(r):return(max(1,np.pi*r**2))#cas discret : une boule de rayon nul, qui contient son centre, est de volume 1.

#Fonction pour le remplissage d'une boule : on suppoe que la norme d'un vecteur de coordonnées est égale à sa norme euclidienne, comme pour les normes 1 et infini


#Fonction de calcul du volume, pour le cas non périodique
##Tient compte seulement de la distance choisie

def vol2D_remp(cen,ray,dim,dist):
 v=0
 for i in range(0,dim[0]):
  for j in range(0,dim[1]):
   if dist([i,j],cen)<=ray:
    v=v+1
 return v

###Dès que delta est non-négatif : l'algorithme ne s'arrête pas en pratique, même pour un petit domaine cf ci-dessus.
###Solution : interrompre la recherche ?
###Dans ce cas, une interruption peut bloquer complètement le remplissage d'une phase, dès que celles qui précèdent auront cessé d'être remplies faute de temps de calcul.









#Nouvelle procédure : toutes les phases sont incluses à chaque tour de boucle, dès que les boules correspondantes satisfont la condition de non-recoupement avec la distance de sécurité delta.
##Un multivolume, comparé aux entrées, décide de la sortie de boucle conjointement au temps d'éxécution.

##On ajoute une entrée : durée d'éxécution pour l'algorithme
##Une autre : périodicité, chaîne de caractères

def RSAA_ph_dist2D(geom,delta,l_ray,par,frac_vol,dim,temps,C_per):
 #volume initial occupé par les phases 1, 2 ... : pour le cas d'une boucle unique ? Pas d'utilité sinon.
 vol_inc= np.zeros(len(frac_vol))
 #sortie : liste unique, dont on peut extraire la liste des centres et rayons concernant une seule phase, en utilisant une liste définie en compréhension.
 liste_ph=[]
 #temps écoulé
 tps=0
 #nombre total de phases
 n_phi=1#len(l_ray)
 #contitions d'arrêt de remplissage : fraction volumique par phase
 C_vol=[True]*n_phi
 #arrêt de la boucle principale :
 Cont_ph=True
 #remplissage
 while Cont_ph:
  for phi in range(0,n_phi):
   l=l_ray[phi]
   #position du centre de la sphère
   cen=[rd.randint(0,dim[0]),rd.randint(0,dim[1])]
   ray=max(0,l(par[phi]))
   #condition : par d'entrecoupement des boules
   C_ent=True
   #test pour toutes les inclusions déjà réalisées : on tient compte des phases délà incluses et on parcourt list_ph
   i=0
   while(C_ent and i<len(liste_ph)):
    if geom[0](liste_ph[i][0],cen)<=delta+liste_ph[i][1]+ray:
     C_ent=False
    #on ajoute la condition qui correspond à la périodicité : boules traversant les six faces du parallélépipède ambiant
    if C_per=='per' :
     if geom[0](liste_ph[i][0]+[dim[0],0],cen)<=delta+liste_ph[i][1]+ray:
      C_ent=False
     if geom[0](liste_ph[i][0]+[-dim[0],0],cen)<=delta+liste_ph[i][1]+ray:
      C_ent=False
     if geom[0](liste_ph[i][0]+[0,dim[1]],cen)<=delta+liste_ph[i][1]+ray:
      C_ent=False
     if geom[0](liste_ph[i][0]+[0,-dim[1]],cen)<=delta+liste_ph[i][1]+ray:
      C_ent=False
    i=i+1
   #On ajoute la nouvelle inclusion, si cela est possible ; on calcule aussi le nouveau volume occupé par les boules
   if C_ent:
    liste_ph.append([cen,ray,phi+1])#décalage entre le numéro de phase et phi, choisi pur parcourir des tableaux
    if C_per=='per':
     vol_inc[phi]=vol_inc[phi]+geom[1](ray)
    else:
     vol_inc[phi]=vol_inc[phi]+vol3D_remp(cen,ray,dim,geom[0])
    ##cas non périodique : volume à calculer avec une autre méthode que geom[1]
   if vol_inc[phi]>frac_vol[phi]*dim[0]*dim[1]:
    C_vol[phi]=False
   tps=tps+1
   #fin de la boucle en phi
  #sortie de la boucle principale
  Cont_ph=any(C_vol) and tps<temps
 return(liste_ph,vol_inc)

##################test

def TestRSAA_ph_dist2D(geom,delta,l_ray,par,frac_vol,dim,temps,C_per):
 #volume initial occupé par les phases 1, 2 ... : pour le cas d'une boucle unique ? Pas d'utilité sinon.
 vol_inc= np.zeros(len(frac_vol))
 #sortie : liste unique, dont on peut extraire la liste des centres et rayons concernant une seule phase, en utilisant une liste définie en compréhension.
 liste_ph=[]
 #temps écoulé
 tps=0
 #nombre total de phases
 n_phi=1#len(l_ray)
 #contitions d'arrêt de remplissage : fraction volumique par phase
 C_vol=[True]*n_phi
 #sortie graphique
 print(n_phi)
 return(liste_ph,vol_inc)


##Remplissage d'un parallélépipède avec les inclusions, d'après une sortie de RSAA_ph_dist2D

def remp2D(ex_rseq,dist,dim,C_per):
 ex=ex_rseq[0]
 A=np.zeros(dim[0]*dim[1])
 A=np.reshape(A,(dim[0],dim[1]))
 #remplissage : attribution d'une phase, point par point du domaine A
 for i in range(0,dim[0]):
  for j in range(0,dim[1]):
   for k in range(0,len(ex)):
    if dist([i,j],ex[k][0])<=ex[k][1]:
     A[i,j]=ex[k][2]
    if C_per=='per':
     if dist([i,j],ex[k][0]+[dim[0],0])<=ex[k][1]:
      A[i,j]=ex[k][2]
     if dist([i,j],ex[k][0]+[-dim[0],0])<=ex[k][1]:
      A[i,j]=ex[k][2]
     if dist([i,j],ex[k][0]+[0,dim[1]])<=ex[k][1]:
      A[i,j]=1+ex[k][2]
     if dist([i,j],ex[k][0]+[0,-dim[1]])<=ex[k][1]:
      A[i,j]=1+ex[k][2]
 return(A)




