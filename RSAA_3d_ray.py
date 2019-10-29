### Paquets a importer ###

import numpy as np
import random as rd


import os


## Remplissage avec des spheres
#Entrees : une norme, une distance minimu entre les inclusions, les lois de probabilite pour les rayons des inbclusions et pour chaque phase, les parametres choisis pour ces lois, une fraction volumique pour chaque phase incluse dans l'espace ambiant initialement vide et un triplet pour le domaine spatial etudie
#Sortie : liste unique, des rayons et centres des inclusions, indexes par le numero de phase.

#Principe : A chaque tour de boucle, on cree une sphere pour une phase donnee, et on verifie que son centre est suffisamment loin de tous les autres centres deja crees.
#Cas periodique ?

##Remarque : l'arret de la procedure principale n'est pas prouve

# A=np.zeros(27)
# A=np.reshape(A,(3,3,3))

#Etape preliminaire : fonctions pour generer des nombres aleatoires dans des boucles. Lois eventuellement tronquees pour obtenir des valeurs positives, on cherche a simuler des rayons

def g_norm(par):
    return(rd.gauss(par[0],par[1]))

#Geometrie choisie avec la premiere variable : formules pour la distance et le volume

##Exemple : geom=[eucl3D,volB] liste de deux fonctions a appeler dans l'algorithme

def eucl3D(a,b):
    return(np.sqrt((a[0]-b[0])**2+(a[1]-b[1])**2+(a[2]-b[2])**2))

def vol3D(r):
    return(max(1,4/3*np.pi*r**3))#cas discret : une boule de rayon nul, qui contient son centre, est de volume 1.

#Fonction de calcul du volume, pour le cas non periodique
##Tient compte seulement de la distance choisie

def vol3D_remp(cen,ray,dim,dist):
    v=0
    for i in range(0,dim[0]):
        for j in range(0,dim[1]):
            for k in range(0,dim[2]):
                if dist([i,j],cen)<=ray:
                    v=v+1
    return v

# #Fonction pour le remplissage d'une boule : on suppoe que la norme d'un vecteur de coordonnees est egale a sa norme euclidienne, comme pour les normes 1 et infini
#
# def remp

# #En pratique : les centres sont distribues selon la loi uniforme de l'espace ambiant.
#
# def RSAA_dist3D(geom,delta,l_ray,par,frac_vol,dim):
#  #volume initial occupe par les phases 1, 2 ... : pour le cas d'une boucle unique ? Pas d'utilite sinon.
#  vol_inc= np.zeros(len(frac_vol))
#  #sortie : liste unique, dont on peut extraire la liste des centres et rayons concernant une seule phase, en utilisant une liste definie en comprehension.
#  liste_ph=[]
#  #nombre total de phases
#  n_phi=len(frac_vol)
#  #remplissage, phase par phase
#  for phi in range(0,n_phi):
#   #sortie de boucle pour le remplissage de la phase courante
#   C_vol=True
#   #Volume initial occupe par la phase i
#   #Vol=0
#   #loi de probabilite parametree suivie par le rayon : appliquer en par[phi]
#   l=l_ray[phi]
#   while C_vol:
#    #position du centre de la sphere
#    cen=[rd.randint(0,dim[0]),rd.randint(0,dim[1]),rd.randint(0,dim[2])]
#    ray=max(0,l(par[phi]))
#    #condition : par d'entrecoupement des boules
#    C_ent=True
#    ##condition supplementaire : boule incluse dans le domaine, cas des normes 1, 2 et infini. Verification faite avant de parcourir liste_ph
#    marge=min(cen[0],cen[1],cen[2],dim[0]-cen[0],dim[1]-cen[1],dim[2]-cen[2])
#    #if marge<ray:
#    # C_ent=False
#    #test pour toutes les inclusions deja realisees : on tient compte des phases dela incluses et on parcourt list_ph
#    i=0
#    while(C_ent and i<len(liste_ph)):
#     if geom[0](liste_ph[i][0],cen)<=delta+liste_ph[i][1]+ray:
#      C_ent=False
#     i=i+1
#    #On ajoute la nouvelle inclusion, si cela est possible
#    if C_ent:
#     liste_ph.append([cen,ray,phi])
#     vol_inc[phi]=vol_inc[phi]+geom[1](ray)
#     print(vol_inc[phi])
#    ##Sortie de boucle : on compare le volume nouvellement occupe par la phase courante a celui prevu en entree
#    if vol_inc[phi]>frac_vol[phi]*dim[0]*dim[1]*dim[2]:
#     C_vol=False
#  return(liste_ph,vol_inc)

# # Exemple d'appel de la fonction RSAA_ray : geometrie euclidienne, 2+1 phases, rayons distribues suivant des lois normales, grille cubique de taille 100
# RSAA_dist3D([eucl3D,vol3D],2,[g_norm,g_norm],[(5,5),(10,1)],[0.15,0.25],[100,100,100])
# ## Exemple miniature
# RSAA_dist3D([eucl3D,vol3D],2,[g_norm,g_norm],[(1,0.5),(2,0.1)],[0.15,0.25],[10,10,10])


###Des que delta est non-negatif : l'algorithme ne s'arrete pas en pratique, meme pour un petit domaine cf ci-dessus.
###Solution : interrompre la recherche ?
###Dans ce cas, une interruption peut bloquer completement le remplissage d'une phase, des que celles qui precedent auront cesse d'etre remplies faute de temps de calcul.

#Nouvelle procedure : toutes les phases sont incluses a chaque tour de boucle, des que les boules correspondantes satisfont la condition de non-recoupement avec la distance de securite delta.
##Un multivolume, compare aux entrees, decide de la sortie de boucle conjointement au temps d'execution.

##On ajoute une entree : duree d'execution pour l'algorithme
##Une autre : periodicite, chaîne de caracteres

def RSAA_ph_dist3D(geom,delta,l_ray,par,frac_vol,dim,temps,C_per):
    # volume initial occupe par les phases 1, 2 ... : pour le cas d'une boucle unique ? Pas d'utilite sinon.
    vol_inc= np.zeros(len(frac_vol))
    # sortie : liste unique, dont on peut extraire la liste des centres et rayons concernant une seule phase, en utilisant une liste definie en comprehension.
    liste_ph=[]
    # temps ecoule
    tps=0
    # nombre total de phases
    n_phi=len(l_ray)
    # contitions d'arret de remplissage : fraction volumique par phase
    C_vol=[True]*n_phi
    # arret de la boucle principale :
    Cont_ph=True
    # remplissage
    while Cont_ph:

        #
        for phi in range(0,n_phi):
            l=l_ray[phi]
            # position du centre de la sphere
            cen=np.array([rd.randint(0,dim[0]),rd.randint(0,dim[1])])
            ray=max(0,l(par[phi]))
            # condition : par d'entrecoupement des boules
            C_ent=True
            # test pour toutes les inclusions deja realisees : on tient compte des phases dela incluses et on parcourt list_ph
            i=0
            while(C_ent and i<len(liste_ph)):
                # une syntaxe avec des tableaux a la place de cen permet d'additionner terme a terme [dim[0],0] et liste_ph[i][0]
                if geom[0](np.array(liste_ph[i][0]),np.array(cen))<=delta+liste_ph[i][1]+ray:
                    C_ent=False
                # on ajoute la condition qui correspond a la periodicite : boules traversant les quatre faces du rectangle ambiant
                if C_per=='per' :
                    if geom[0](liste_ph[i][0]+[dim[0],0,0],cen)<=delta+liste_ph[i][1]+ray:
                        C_ent=False
                    if geom[0](liste_ph[i][0]+[-dim[0],0,0],cen)<=delta+liste_ph[i][1]+ray:
                        C_ent=False
                    if geom[0](liste_ph[i][0]+[0,dim[1],0],cen)<=delta+liste_ph[i][1]+ray:
                        C_ent=False
                    if geom[0](liste_ph[i][0]+[0,-dim[1],0],cen)<=delta+liste_ph[i][1]+ray:
                        C_ent=False
                    if geom[0](liste_ph[i][0]+[0,0,dim[2]],cen)<=delta+liste_ph[i][1]+ray:
                        C_ent=False
                    if geom[0](liste_ph[i][0]+[0,0-dim[2]],cen)<=delta+liste_ph[i][1]+ray:
                        C_ent=False
                i=i+1

            # On ajoute la nouvelle inclusion, si cela est possible ; on calcule aussi le nouveau volume occupe par les boules
            if C_ent:
                # ajout de l'inclusion, decalage entre le numero de phase et phi, choisi pur parcourir des tableaux
                liste_ph.append([cen,ray,phi+1])
                # calcul du nouveau volume occupe par la pĥase phi
                if C_per=='per':
                    vol_inc[phi]=vol_inc[phi]+geom[1](ray)
                else:
                    ## cas non periodique : volume a calculer avec une autre methode que geom[1]
                    vol_inc[phi]=vol_inc[phi]+vol3D_remp(cen,ray,dim,geom[0])
            # continuation du remplissage : condition sur la fraction volumique
            if vol_inc[phi]>frac_vol[phi]*dim[0]*dim[1]*dim[2]:
                C_vol[phi]=False

            tps=tps+1
            # fin de la boucle en phi

        # sortie de la boucle principale
        Cont_ph=any(C_vol) and tps<temps

    return(liste_ph,vol_inc)


# RSAA_ph_dist3D([eucl3D,vol3D],2,[g_norm,g_norm],[(5,5),(10,1)],[0.15,0.25],[100,100,100],10000,False)
#
# RSAA_ph_dist3D([eucl3D,vol3D],2,[g_norm,g_norm],[(5,5),(10,1)],[0.15,0.25],[100,100,100],10000,True)
#
# #RSAA_ph_dist3D([eucl3D,vol3D],2,[g_norm,g_norm],[(5,5),(10,1)],[0.15,0.25],[100,100,100],100000,True)
# ##... environ quatre minutes, fractions volumiques array([ 161719.35151185,  138658.97646631])), peu de boules "1" et beaucoup de boules "0" de rayon nul
#
# RSAA_ph_dist3D([eucl3D,vol3D],2,[g_norm,g_norm],[(5,2),(7,1)],[0.15,0.25],[100,100,100],100000,'per')
# #environ deux minutes, array([ 141891.86874144,  113682.71909738])),
#
# ##Exemple miniature
#
# RSAA_ph_dist3D([eucl3D,vol3D],2,[g_norm,g_norm],[(1,0.5),(2,0.1)],[0.15,0.25],[10,10,10],100,'per')
#
# RSAA_ph_dist3D([eucl3D,vol3D],4,[g_norm,g_norm,g_norm],[(3,0.5),(8,0.2),(2,0.1)],[0.15,0.25,0.1],[30,50,10],10000,'per')
