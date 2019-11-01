### Encodage ###


### Paquets a importer ###

#if __name__=='__main__':
import numpy as np
import random as rd

#import pylab as pl
#from pylab import *


#import os


## Remplissage avec des spheres
#Entrees : une norme, une distance minimu entre les inclusions, les lois de probabilite pour les rayons des inclusions et pour chaque phase, les parametres choisis pour ces lois, une fraction volumique pour chaque phase incluse dans l'espace ambiant initialement vide et un triplet pour le domaine spatial etudie
#Sortie : liste unique, des rayons et centres des inclusions, indexes par le numero de phase.

#Principe : A chaque tour de boucle, on cree une sphere pour une phase donnee, et on verifie que son centre est suffisamment loin de tous les autres centres deja crees.
#Cas periodique ?

##Remarque : l'arret de la procedure principale n'est pas prouve

#Etape preliminaire : fonctions pour generer des nombres aleatoires dans des boucles. Lois eventuellement tronquees pour obtenir des valeurs positives, on cherche a simuler des rayons

def g_norm(par):

    return(rd.gauss(par[0],par[1]))

def unif(par):

    return(par[0] + rd.random()*(par[1] - par[0]))

#Geometrie choisie avec la premiere variable : formules pour la distance et le volume

##Exemple : geom=[eucl2D,volB] liste de deux fonctions a appeler dans l'algorithme

def eucl2D(a,b):
    return(np.sqrt((a[0]-b[0])**2+(a[1]-b[1])**2))

def vol2D(r):
    return(max(1,np.pi*r**2))#cas discret : une boule de rayon nul, qui contient son centre, est de volume 1.

def vol2D_size(r):
    return(np.pi*r**2)#cas discret : une boule de rayon nul, qui contient son centre, est de volume 1.


#Fonction pour le remplissage d'une boule : on suppoe que la norme d'un vecteur de coordonnees est egale a sa norme euclidienne, comme pour les normes 1 et infini


#Fonction de calcul du volume, pour le cas non periodique
##Tient compte seulement de la distance choisie

def vol2D_remp(cen,ray,dim,dist):
    v=0
    for i in range(0,dim[0]):
            for j in range(0,dim[1]):
                    if dist([i,j],cen)<=ray:
                        v=v+1
    return v

###Des que delta est non-negatif : l'algorithme ne s'arrete pas en pratique, meme pour un petit domaine cf ci-dessus.
###Solution : interrompre la recherche ?
###Dans ce cas, une interruption peut bloquer completement le remplissage d'une phase, des que celles qui precedent auront cesse d'etre remplies faute de temps de calcul.









#Nouvelle procedure : toutes les phases sont incluses a chaque tour de boucle, des que les boules correspondantes satisfont la condition de non-recoupement avec la distance de securite delta.
##Un multivolume, compare aux entrees, decide de la sortie de boucle conjointement au temps d'execution.

##On ajoute une entree : duree d'execution pour l'algorithme
##Une autre : periodicite, chaîne de caracteres

def RSAA_ph_dist2D(geom,delta,l_ray,par,frac_vol,dim,temps,C_per):
    #volume initial occupe par les phases 1, 2 ... : pour le cas d'une boucle unique ? Pas d'utilite sinon.
    vol_inc= np.zeros(len(frac_vol))
    #sortie : liste unique, dont on peut extraire la liste des centres et rayons concernant une seule phase, en utilisant une liste definie en comprehension.
    liste_ph=[]
    #temps ecoule
    tps=0
    #nombre total de phases
    n_phi=len(l_ray)
    #contitions d'arret de remplissage : fraction volumique par phase
    C_vol=[True]*n_phi
    #arret de la boucle principale :
    Cont_ph=True
    #remplissage
    while Cont_ph:
        for phi in range(0,n_phi):
            l=l_ray[phi]
            #position du centre de la sphere
            cen=np.array([rd.randint(0,dim[0]),rd.randint(0,dim[1])])
            ray=max(0,l(par[phi]))
            #condition : par d'entrecoupement des boules
            C_ent=True
            #test pour toutes les inclusions deja realisees : on tient compte des phases dela incluses et on parcourt list_ph
            i=0
            while(C_ent and i<len(liste_ph)):
                if geom[0](np.array(liste_ph[i][0]),np.array(cen))<=delta+liste_ph[i][1]+ray:#une syntaxe avec des tableaux a la place de cen permet d'additionner terme a terme [dim[0],0] et liste_ph[i][0]
                    C_ent=False
                #on ajoute la condition qui correspond a la periodicite : boules traversant les quatre faces du rectangle ambiant
                if C_per=='per' :
                    if geom[0](np.array(liste_ph[i][0])+np.array([dim[0],0]),cen)<=delta+liste_ph[i][1]+ray:
                        C_ent=False
                    if geom[0](np.array(liste_ph[i][0])+np.array([-dim[0],0]),cen)<=delta+liste_ph[i][1]+ray:
                        C_ent=False
                    if geom[0](np.array(liste_ph[i][0])+np.array([0,dim[1]]),cen)<=delta+liste_ph[i][1]+ray:
                        C_ent=False
                    if geom[0](np.array(liste_ph[i][0])+np.array([0,-dim[1]]),cen)<=delta+liste_ph[i][1]+ray:
                        C_ent=False
                i=i+1
            #On ajoute la nouvelle inclusion, si cela est possible ; on calcule aussi le nouveau volume occupe par les boules
            if C_ent:
                liste_ph.append([cen,ray,phi+1])#decalage entre le numero de phase et phi, choisi pur parcourir des tableaux
                if C_per=='per':
                    vol_inc[phi]=vol_inc[phi]+geom[1](ray)
                else:
                    vol_inc[phi]=vol_inc[phi]+vol2D_remp(cen,ray,dim,geom[0])
                ##cas non periodique : volume a calculer avec une autre methode que geom[1]
            if vol_inc[phi]>frac_vol[phi]*dim[0]*dim[1]:
                C_vol[phi]=False
            tps=tps+1
            #fin de la boucle en phi
        #sortie de la boucle principale
        Cont_ph=any(C_vol) and tps<temps
    return(liste_ph,vol_inc)

# RSAA dans une cellule elementaire exploitable avec rom_diffeo_dhom

# Une autre fonction, adaptee au choix d'une cellule elementaire

def RSAA_ph_eucl_cell(delta, l_ray, par, frac_vol, xyinf, size, temps):

    # volume initial occupe par les phases 1, 2 ... : pour le cas d'une boucle unique ? Pas d'utilite sinon.
    vol_inc= np.zeros(len(frac_vol))

    # bornes de la cellule
    xinf = xyinf[0]
    yinf = xyinf[1]
    #
    xsup = xinf + size
    ysup = yinf + size

    # sortie : liste unique, dont on peut extraire la liste des centres et rayons concernant une seule phase, en utilisant une liste definie en comprehension.
    liste_ph=[]

    # nombre total de phases
    n_phi=len(l_ray)

    # une iteration : placement d'au plus une inclusion par phase
    for tps in range(temps):

        # pour chaque phase solide : on tente de placer une inclusion
        for phi in range(0,n_phi):

            # condition pour placer l'inclusion : fraction volumique imposee
            if vol_inc[phi] < frac_vol[phi]*size**2:

                # loi  de probabilite pour le rayon a venir et tirage de ce rayon
                l = l_ray[phi]
                # ray = max(0,l(par[phi]))
                ray = abs(l(par[phi]))

                # position du centre de la sphere tiree uniformement
                cen = np.array([rd.uniform(xinf, xsup), rd.uniform(yinf, ysup)])

                # condition : par d'entrecoupement des boules
                C_ent=True
                # test pour toutes les inclusions deja realisees : on tient compte des phases dela incluses et on parcourt list_ph
                i = 0
                # test sur l'entrecoupement
                while(C_ent and i<len(liste_ph)):

                    # une syntaxe avec des tableaux a la place de cen permet d'additionner terme a terme [dim[0],0] et liste_ph[i][0]
                    cen_test = liste_ph[i][0]
                    ray_test = liste_ph[i][1]


                    if eucl2D(cen_test, cen) <= delta + ray_test + ray:
                        C_ent=False

                    # # on ajoute la condition qui correspond a la periodicite : boules traversant les quatre faces du rectangle ambiant
                    # if eucl2D(cen_test+[size,0], cen) <= delta + ray_test + ray:
                    #     C_ent=False
                    # if eucl2D(cen_test+[-size,0], cen) <= delta + ray_test + ray:
                    #     C_ent=False
                    # if eucl2D(cen_test+[0,size], cen) <= delta + ray_test + ray:
                    #     C_ent=False
                    # if eucl2D(cen_test+[0,-size], cen) <= delta + ray_test + ray:
                    #     C_ent=False
                    # ## --- fonctionne --- ##

                    # on ajoute la condition qui correspond a la periodicite : boules traversant les quatre faces du rectangle ambiant
                    for a in range(2):
                        for h in [-size, size]:

                            vect = np.zeros(2)
                            vect[a] = h

                            cen_test_per = cen_test + vect
                            # cen_test_per[a] = cen_test[a] + vect

                            if eucl2D(cen_test_per, cen) <= delta + ray_test + ray:
                                C_ent=False
                    ## --- fonctionne : plus de chevauchement --- ##

                    # on passe a l'inclusion suivante pour testerle recoupement
                    i = i + 1

                # On ajoute la nouvelle inclusion, si cela est possible ; on calcule aussi le nouveau volume occupe par les boules
                if C_ent:
                    # ajout de l'inclusion, decalage entre le numero de phase et phi, choisi pur parcourir des tableaux
                    liste_ph.append([cen, ray, phi+1])
                    # calcul du nouveau volume occupe par la pĥase phi
                    vol_inc[phi] = vol_inc[phi] + vol2D_size(ray)

        # fraction volumique occupee par les phases solides
        total_frac_vol = vol_inc/size**2

    return(liste_ph, total_frac_vol)


# Une fonction pour ellipses : preliminaires

## ellipse : [[cen], [gax, dil], theta]

def dist2D_ellell(cen_rdil_th_1, cen_rdil_th_2):

    cen_1 = cen_rdil_th_1[0]
    cen_2 = cen_rdil_th_2[0]

    dist2D = np.sqrt((cen_1[0] - cen_2[0])**2 + (cen_1[1] - cen_2[1])**2)

    return(dist2D)

def vol2D_ellell(cen_rdil_th):

    gax = cen_rdil_th[1][0]
    dil = cen_rdil_th[1][1]

    vol2D = np.pi*dil*(gax**2)

    return(vol2D)

# Une fonction pour ellipses : RSAA avec theta et pax/gax tires selon une loi uniforme

def RSAA_ph_ell_cell(delta, l_ray, par, frac_vol, xyinf, size, temps):

    # volume initial occupe par les phases 1, 2 ... : pour le cas d'une boucle unique ? Pas d'utilite sinon.
    vol_inc= np.zeros(len(frac_vol))

    # bornes de la cellule
    xinf = xyinf[0]
    yinf = xyinf[1]
    #
    xsup = xinf + size
    ysup = yinf + size

    # sortie : liste unique, dont on peut extraire la liste des centres et rayons concernant une seule phase, en utilisant une liste definie en comprehension.
    liste_ph=[]

    # nombre total de phases
    n_phi=len(l_ray)

    # une iteration : placement d'au plus une inclusion par phase
    for tps in range(temps):

        # pour chaque phase solide : on tente de placer une inclusion
        for phi in range(0,n_phi):

            # condition pour placer l'inclusion : fraction volumique imposee
            if vol_inc[phi] < frac_vol[phi]*size**2:

                # loi  de probabilite pour le rayon a venir et tirage de ce rayon
                l = l_ray[phi]
                gax = abs(l(par[phi]))

                # position du centre de la sphere tiree uniformement
                cen = np.array([rd.uniform(xinf, xsup), rd.uniform(yinf, ysup)])
                # rapport du petit axe au grand axe
                dil = rd.uniform(0.5, 1)
                # angle d'inclinaison du grand axe par rapport a l'axe des abscisses
                theta = rd.uniform(0, 2*np.pi)

                # ellipse courante
                ell = [cen, [gax, dil], theta]

                # condition : par d'entrecoupement des boules
                C_ent=True
                # test pour toutes les inclusions deja realisees : on tient compte des phases dela incluses et on parcourt list_ph
                i = 0
                # test sur l'entrecoupement
                while(C_ent and i<len(liste_ph)):

                    # inclusion a tester
                    ell_test = liste_ph[i][0]
                    # composantes a extraire, pour translater vers les cellules voisines et evaluer des distances
                    cen_test = ell_test[0]
                    gax_test = ell_test[1][0]
                    dil_test = ell_test[1][1]
                    theta_test = ell_test[2]

                    # distance entre l'inclusion courante et l'inclusion testee
                    if dist2D_ellell(ell_test, ell) <= delta + gax_test + gax:
                        C_ent=False

                    # on ajoute la condition qui correspond a la periodicite : boules traversant les quatre faces du rectangle ambiant
                    for a in range(2):
                        for h in [-size, size]:

                            vect = np.zeros(2)
                            vect[a] = h

                            # on deplace le centre de l'ellipse testee
                            cen_test_per = cen_test + vect
                            ell_test_per = [cen_test_per, [gax_test, dil_test], theta_test]

                            if dist2D_ellell(ell_test_per, ell) <= delta + gax_test + gax:
                                C_ent=False

                    # on passe a l'inclusion suivante pour testerle recoupement
                    i = i + 1


                # On ajoute la nouvelle inclusion, si cela est possible ; on calcule aussi le nouveau volume occupe par les boules
                if C_ent:
                    # ajout de l'inclusion, decalage entre le numero de phase et phi, choisi pur parcourir des tableaux
                    liste_ph.append([ell, phi+1])
                    # if tps <= 5:
                    #     print('%'*60)
                    #     print('ellipse courante :', ell)
                    #     print('-'*60)
                    #     print('derniere inclusion :', liste_ph[len(liste_ph)-1])
                    # calcul du nouveau volume occupe par la pĥase phi
                    vol_inc[phi] = vol_inc[phi] + vol2D_ellell(ell)

                    # if tps <= 5:
                    #     print('-'*60)
                    #     print('premiere inclusion apres vol2D :', liste_ph[0])
                    #     print('%'*60)

        # fraction volumique occupee par les phases solides
        total_frac_vol = vol_inc/size**2

    return(liste_ph, total_frac_vol)


















##Remplissage d'un rectangle avec les inclusions, d'apres une sortie de RSAA_ph_dist2D

def remp2D(ex_rseq,dist,dim,C_per):
    ex=ex_rseq[0]
    #z=0
    A=np.zeros(dim[0]*dim[1])
    A=np.reshape(A,(dim[0],dim[1]))
    #remplissage : attribution d'une phase, point par point du domaine A
    for i in range(0,dim[0]):
        for j in range(0,dim[1]):
            for k in range(0,len(ex)):
                if dist([i,j],ex[k][0])<=ex[k][1]:
                    A[i,j]=ex[k][2]
                if C_per=='per':
                    #z=z+1
                    #print(dist([i,j],ex[k][0]+[dim[0],0]),ex[k][1])
                    if dist(np.array([i,j]),ex[k][0]+np.array([dim[0],0]))<=ex[k][1]:
                        A[i,j]=ex[k][2]
                    if dist(np.array([i,j]),ex[k][0]+np.array([-dim[0],0]))<=ex[k][1]:
                        A[i,j]=ex[k][2]
                    if dist(np.array([i,j]),ex[k][0]+np.array([0,dim[1]]))<=ex[k][1]:
                        A[i,j]=ex[k][2]
                    if dist(np.array([i,j]),ex[k][0]+np.array([0,-dim[1]]))<=ex[k][1]:
                        A[i,j]=ex[k][2]
    return(A)#[A,z])

#Optimisation : une fonction vectorisee, peu d'interet pour de grandes dimensions

#Une fonction intermediaire : phase au pixel I [i,j]

def phase_pt(I,liste_ph,n_phi,dist,dim,C_per):
    phase=0
    C_ont=True
    phi=1
    i=I[0]
    j=I[1]
    while C_ont and phi<=n_phi:
        L=[x for x in liste_ph if x[2]==phi]
        l=len(L)
        A=np.array(L)
        #on extrait les vecteurs correspondant aux rayons et aux centres
        cen=A[:,0]
        ray=A[:,1]
        #fonction distance vectorisee
        def dist_pt(pt):
            return(dist([i,j],pt))
        v_dist_ij_pt=np.vectorize(dist_pt)
        #calcul des distances aux centres :
        dist_cen=v_dist_ij_pt(cen)
        if any(dist_cen<=ray):
            phase=phi
            C_ont=False
        #cas periodique
        if C_per=='per':
            ##------------------------- on translate les boules de [dim[],] etc puis on verifie l'appartenance du point courant I a leurs images
            if C_ont:
                B=A+[[dim[0],0],0,0]
                cen=B[:,0]
                ray=B[:,1]
                #print(cen)
                #print(ray)
                dist_cen=v_dist_ij_pt(cen)
                #print(dist_cen)
            if any(dist_cen<=ray):
                phase=phi
                C_ont=False
            if C_ont:
                B=A+[[-dim[0],0],0,0]
                cen=B[:,0]
                ray=B[:,1]
                dist_cen=v_dist_ij_pt(cen)
            if any(dist_cen<=ray):
                phase=phi
                C_ont=False
            if C_ont:
                B=A+[[0,dim[1]],0,0]
                cen=B[:,0]
                ray=B[:,1]
                dist_cen=v_dist_ij_pt(cen)
            if any(dist_cen<=ray):
                phase=phi
                C_ont=False
            if C_ont:
                B=A+[[0,-dim[1]],0,0]
                cen=B[:,0]
                ray=B[:,1]
                dist_cen=v_dist_ij_pt(cen)
            if any(dist_cen<=ray):
                phase=phi
                C_ont=False
        ##------------------------- fin pour la periodicite ----------------------------##
        phi=phi+1
    return(phase)

#Reformulation : indice k unique pour le parcours de la matrice spatiale dans l'ordre lexicographique, eventuellement a lancer sur des treads independants

def parc_liste(k_lex,L,n_phi,dist,dim,C_per):
    i=k_lex//dim[1]
    j=k_lex-i*dim[1]
    a_ph=phase_pt([i,j],L,n_phi,dist,dim,C_per)
    return(a_ph)

#Remplissage, boucle avec des instructions vectorisees, un seul thread

def Vremp2D(ex_rseq,dist,dim,C_per):
    A=np.arange(dim[0]*dim[1])
    L=ex_rseq[0]#liste des centre-rayon-phase ; on oublie les saturations respectives des phases
    M=np.array(L)
    n_phi=max(M[:,2])
    #boucle principale
    for k in range(0,len(A)):
        A[k]=parc_liste(k,L,n_phi,dist,dim,C_per)
        #sortie sous la forme d'une matrice : pixels
    A=np.reshape(A,(dim[0],dim[1]))
    return(A)



#version parallelisee


#if __name__=='__main__':
import multiprocessing
from joblib import Parallel, delayed


def ParVremp2D(ex_rseq,dist,dim,C_per):
    A=np.arange(dim[0]*dim[1])
    L=ex_rseq[0]#liste des centre-rayon-phase ; on oublie les saturations respectives des phases
    M=np.array(L)
    n_phi=max(M[:,2])
    #definition et execution des threads : huit, pour MECALAC
    ##
    num_cores=multiprocessing.cpu_count()
    ##
    #boucle principale, parallelisee
    ##
    A=Parallel(n_jobs=num_cores)(delayed(parc_liste)(k,L,n_phi,dist,dim,C_per)for k in range(0,len(A)))
    ##
    A=np.reshape(A,(dim[0],dim[1]))
    return(A)
