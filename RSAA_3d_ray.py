### Paquets a importer ###

import numpy as np
import random as rd


import os



# # caracteristiques de la cellule
#
# xinf = -1.
# yinf = -1.
# zinf = -1.
#
# size = 2.
#
# xsup = xinf + size
# ysup = yinf + size
# zsup = zinf + size

## Remplissage avec des spheres


# Etape preliminaire : fonctions pour generer des nombres aleatoires dans des boucles. Lois eventuellement tronquees pour obtenir des valeurs positives, on cherche a simuler des rayons

def g_norm(par):
    return(rd.gauss(par[0],par[1]))

def unif(par):
    # return(rd.uniform(par[0],par[1]))

    return(par[0] + rd.random()*(par[1] - par[0]))

# Geometrie choisie avec la premiere variable : formules pour la distance et le volume

## Exemple : geom=[eucl3D,volB] liste de deux fonctions a appeler dans l'algorithme

def eucl3D(a,b):
    return(np.sqrt((a[0]-b[0])**2+(a[1]-b[1])**2+(a[2]-b[2])**2))

def vol3D(r):
    return(4/3*np.pi*r**3)#cas discret : une boule de rayon nul, qui contient son centre, est de volume 1.

# Fonction de calcul du volume, pour le cas non periodique
## Tient compte seulement de la distance choisie

def vol3D_remp(cen,ray,dim,dist):
    v=0
    for i in range(0,dim[0]):
        for j in range(0,dim[1]):
            for k in range(0,dim[2]):
                if dist([i,j],cen)<=ray:
                    v=v+1
    return v




# Nouvelle procedure : toutes les phases sont incluses a chaque tour de boucle, des que les boules correspondantes satisfont la condition de non-recoupement avec la distance de securite delta.
## Un multivolume, compare aux entrees, decide de la sortie de boucle conjointement au temps d'execution.


# Une autre fonction, adaptee au choix d'une cellule elementaire

def RSAA_ph_eucl_cell(delta, l_ray, par, frac_vol, xyzinf, size, temps, C_per):

    # volume initial occupe par les phases 1, 2 ... : pour le cas d'une boucle unique ? Pas d'utilite sinon.
    vol_inc= np.zeros(len(frac_vol))

    # bornes de la cellule
    xinf = xyzinf[0]
    yinf = xyzinf[1]
    zinf = xyzinf[2]
    #
    xsup = xinf + size
    ysup = yinf + size
    zsup = zinf + size

    # sortie : liste unique, dont on peut extraire la liste des centres et rayons concernant une seule phase, en utilisant une liste definie en comprehension.
    liste_ph=[]

    # nombre total de phases
    n_phi=len(l_ray)

    # une iteration : placement d'au plus une inclusion par phase
    for tps in range(temps):

        # pour chaque phase solide : on tente de placer une inclusion
        for phi in range(0,n_phi):

            # condition pour placer l'inclusion : fraction volumique imposee
            if vol_inc[phi] < frac_vol[phi]*size**3:

                # loi  de probabilite pour le rayon a venir et tirage de ce rayon
                l = l_ray[phi]
                # ray = max(0,l(par[phi]))
                ray = abs(l(par[phi]))

                # position du centre de la sphere tiree uniformement
                cen = np.array([rd.uniform(xinf, xsup), rd.uniform(yinf, ysup), rd.uniform(zinf, zsup)])

                # condition : par d'entrecoupement des boules
                C_ent=True
                # test pour toutes les inclusions deja realisees : on tient compte des phases dela incluses et on parcourt list_ph
                i = 0
                # test sur l'entrecoupement
                while(C_ent and i<len(liste_ph)):
                    # une syntaxe avec des tableaux a la place de cen permet d'additionner terme a terme [dim[0],0] et liste_ph[i][0]
                    if eucl3D(liste_ph[i][0],cen)<=delta+liste_ph[i][1]+ray:
                        C_ent=False
                    # on ajoute la condition qui correspond a la periodicite : boules traversant les quatre faces du rectangle ambiant
                    if C_per=='per' :
                        if eucl3D(liste_ph[i][0]+[size,0,0],cen)<=delta+liste_ph[i][1]+ray:
                            C_ent=False
                        if eucl3D(liste_ph[i][0]+[-size,0,0],cen)<=delta+liste_ph[i][1]+ray:
                            C_ent=False
                        if eucl3D(liste_ph[i][0]+[0,size,0],cen)<=delta+liste_ph[i][1]+ray:
                            C_ent=False
                        if eucl3D(liste_ph[i][0]+[0,-size,0],cen)<=delta+liste_ph[i][1]+ray:
                            C_ent=False
                        if eucl3D(liste_ph[i][0]+[0,0,size],cen)<=delta+liste_ph[i][1]+ray:
                            C_ent=False
                        if eucl3D(liste_ph[i][0]+[0,0,-size],cen)<=delta+liste_ph[i][1]+ray:
                            C_ent=False
                    i = i+1

                # On ajoute la nouvelle inclusion, si cela est possible ; on calcule aussi le nouveau volume occupe par les boules
                if C_ent:
                    # ajout de l'inclusion, decalage entre le numero de phase et phi, choisi pur parcourir des tableaux
                    liste_ph.append([cen, ray, phi+1])
                    # calcul du nouveau volume occupe par la pĥase phi
                    if C_per == 'per':
                        vol_inc[phi] = vol_inc[phi] + vol3D(ray)
                    else:
                        ## cas non periodique : volume a calculer avec une autre methode que geom[1]
                        vol_inc[phi] = vol_inc[phi] + vol3D_remp(cen, ray, dim, eucl3D)

        # fraction volumique occupee par les phases solides
        total_frac_vol = vol_inc/size**3

    return(liste_ph, total_frac_vol)









# Une fonction pour spheres, et cylindres dans les trois directions orthogonales d'espace ; periodicite supposee, geometrie euclidienne

def dist3D_cylcyl(ax_dir_1, ax_dir_2):

    dir_1 = ax_dir_1[1]
    dir_2 = ax_dir_2[1]

    ax_1 = ax_dir_1[0]
    ax_2 = ax_dir_2[0]

    if dir_1 == dir_2:
        # on annule la coordonnee variable
        ax_1[dir_1] = 0
        ax_2[dir_2] = 0
        # on en deduit la distance entre les deux axes paralleles
        dist3D = np.sqrt((ax_1[0]-ax_2[0])**2 + (ax_1[1]-ax_2[1])**2 + (ax_1[2]-ax_2[2])**2)

    else :
        # on amene les points courants des axes sur la perpendiculaire commune aux axes
        ax_1[dir_1] = ax_2[dir_1]
        ax_2[dir_2] = ax_1[dir_2]
        # on en deduit la distance entre les deux axes paralleles
        dist3D = np.sqrt((ax_1[0]-ax_2[0])**2 + (ax_1[1]-ax_2[1])**2 + (ax_1[2]-ax_2[2])**2)

    return(dist3D)


def dist3D_cylsph(ax_dir, cen):

    ax = ax_dir[0]
    dir = ax_dir[1]

    # projete orthogonal du centre de la spĥhere sur l'axe du cylindre
    proj = ax
    proj[dir] = cen[dir]

    # distance du centre de la sphere a l'axe du cylindre
    dist3D = np.sqrt((cen[0]-proj[0])**2 + (cen[1]-proj[1])**2 + (cen[2]-proj[2])**2)

    return(dist3D)




def RSAA_ph_cylsph_cell(delta, l_ray, par, frac_vol, xyzinf, size, temps):

    # volume initial occupe par les phases 1, 2 ... : pour le cas d'une boucle unique ? Pas d'utilite sinon.
    vol_inc= np.zeros(len(frac_vol))

    # bornes de la cellule
    xinf = xyzinf[0]
    yinf = xyzinf[1]
    zinf = xyzinf[2]
    #
    xsup = xinf + size
    ysup = yinf + size
    zsup = zinf + size

    # sortie : liste unique, dont on peut extraire la liste des centres et rayons concernant une seule phase, en utilisant une liste definie en comprehension.
    liste_ph=[]

    # nombre total de phases
    n_phi = 2

    # une iteration : placement d'au plus une inclusion par phase
    for tps in range(temps):

        # ---------------------- ajout d'un cylindre : phase solide 1 ---------------------- #
        phi = 0
        # condition pour placer l'inclusion : fraction volumique imposee
        if vol_inc[phi] < frac_vol[phi]*size**3:

            # loi  de probabilite pour le rayon a venir et tirage de ce rayon
            l = l_ray[phi]
            # ray = max(0,l(par[phi]))
            ray = abs(l(par[phi]))

            # point arbitraire de l'axe, coordonnees a exploiter
            ax_pt = np.array([rd.uniform(xinf, xsup), rd.uniform(yinf, ysup), rd.uniform(zinf, zsup)])
            # tirage d'une direction
            dir = rd.randint(0, 2)
            # caracteristiques de l'axes
            ax_dir = [ax_pt, dir]

            # condition : par d'entrecoupement des boules
            C_ent=True
            # test pour toutes les inclusions deja realisees : on tient compte des phases dela incluses et on parcourt list_ph
            i = 0
            # test sur l'entrecoupement
            while(C_ent and i<len(liste_ph)):
                # centre de sphere ou axe de cylindre dirige
                cen_or_axdir = liste_ph[i][0]
                ray_test = liste_ph[i][1]
                phase = liste_ph[i][2]
                # rencontre avec une sphere et ses six translatees
                if phase == 2:
                    # centre de la sphere a tester
                    cen_test = cen_or_axdir
                    # une syntaxe avec des tableaux a la place de cen permet d'additionner terme a terme [dim[0],0] et liste_ph[i][0]
                    if dist3D_cylsph(ax_dir, cen_test) <= delta + ray_test + ray:
                        C_ent=False
                    # on ajoute la condition qui correspond a la periodicite : boules traversant les quatre faces du rectangle ambiant
                    if dist3D_cylsph(ax_dir, cen_test + [size, 0, 0]) <= delta + ray_test + ray:
                        C_ent=False
                    if dist3D_cylsph(ax_dir, cen_test + [-size, 0, 0]) <= delta + ray_test + ray:
                        C_ent=False
                    if dist3D_cylsph(ax_dir, cen_test + [0, size, 0]) <= delta + ray_test + ray:
                        C_ent=False
                    if dist3D_cylsph(ax_dir, cen_test + [0, -size, 0]) <= delta + ray_test + ray:
                        C_ent=False
                    if dist3D_cylsph(ax_dir, cen_test + [0, 0, size]) <= delta + ray_test + ray:
                        C_ent=False
                    if dist3D_cylsph(ax_dir, cen_test + [0, 0, -size]) <= delta + ray_test + ray:
                        C_ent=False

                else:
                    ax_test = cen_or_axdir[0]
                    dir_test = cen_or_axdir[1]
                    # test sans periodicite
                    if dist3D_cylcyl(ax_dir, cen_or_axdir) <= delta + ray_test + ray:
                        C_ent=False
                    # periodicite :
                    for a in range(len(ax_test)):
                        # quatre translations du cylindre rencontre
                        if a != dir_test:
                            for vect in [-size, size]:
                                ax_test_per = ax_test
                                ax_test_per[a] += vect
                                axdir_test_per = [ax_test_per, dir_test]
                                # test apres translation du cylindre rencontre
                                if dist3D_cylcyl(ax_dir, axdir_test_per) <= delta + ray_test + ray:
                                    C_ent=False

                # on continue le parcours de la liste des inclusions existantes
                i = i + 1

            # On ajoute la nouvelle inclusion, si cela est possible ; on calcule aussi le nouveau volume occupe par les boules
            if C_ent:
                # ajout de l'inclusion, decalage entre le numero de phase et phi, choisi pur parcourir des tableaux
                liste_ph.append([ax_dir, ray, phi+1])
                # calcul du nouveau volume occupe par la pĥase phi
                vol_inc[phi] = vol_inc[phi] + size*np.pi*ray**2

        # ---------------------- ajout d'une sphere : phase solide 2 ---------------------- #
        phi = 1
        # condition pour placer l'inclusion : fraction volumique imposee
        if vol_inc[phi] < frac_vol[phi]*size**3:

            # loi  de probabilite pour le rayon a venir et tirage de ce rayon
            l = l_ray[phi]
            # ray = max(0,l(par[phi]))
            ray = abs(l(par[phi]))

            # position du centre de la sphere tiree uniformement
            cen = np.array([rd.uniform(xinf, xsup), rd.uniform(yinf, ysup), rd.uniform(zinf, zsup)])

            # condition : par d'entrecoupement des boules
            C_ent=True
            # test pour toutes les inclusions deja realisees : on tient compte des phases dela incluses et on parcourt list_ph
            i = 0
            # test sur l'entrecoupement
            while(C_ent and i<len(liste_ph)):

                # centre de sphere ou axe de cylindre dirige
                cen_or_axdir = liste_ph[i][0]
                ray_test = liste_ph[i][1]
                phase = liste_ph[i][2]
                # rencontre avec une sphere et ses six translatees
                if phase == 2:
                    # une syntaxe avec des tableaux a la place de cen permet d'additionner terme a terme [dim[0],0] et liste_ph[i][0]
                    if eucl3D(cen_or_axdir, cen)<=delta + liste_ph[i][1] + ray:
                        C_ent=False
                    # on ajoute la condition qui correspond a la periodicite : boules traversant les quatre faces du rectangle ambiant
                    if eucl3D(cen_or_axdir + [size,0,0], cen)<=delta + ray_test + ray:
                        C_ent=False
                    if eucl3D(cen_or_axdir + [-size,0,0], cen)<=delta + ray_test + ray:
                        C_ent=False
                    if eucl3D(cen_or_axdir + [0,size,0], cen)<=delta + ray_test + ray:
                        C_ent=False
                    if eucl3D(cen_or_axdir + [0,-size,0], cen)<=delta + ray_test + ray:
                        C_ent=False
                    if eucl3D(cen_or_axdir + [0,0,size], cen)<=delta + ray_test + ray:
                        C_ent=False
                    if eucl3D(cen_or_axdir + [0,0,-size], cen)<=delta + ray_test + ray:
                        C_ent=False

                else :
                    ax_test = cen_or_axdir[0]
                    dir_test = cen_or_axdir[1]
                    # test sans periodicite
                    if dist3D_cylsph(cen_or_axdir, cen) <= delta + ray_test + ray:
                        C_ent=False
                    # periodicite :
                    for a in range(len(ax_test)):
                        # quatre translations du cylindre rencontre
                        if a != dir_test:
                            for vect in [-size, size]:
                                ax_test_per = ax_test
                                ax_test_per[a] += vect
                                axdir_test_per = [ax_test_per, dir_test]
                                # test apres translation du cylindre rencontre
                                if dist3D_cylsph(axdir_test_per, cen) <= delta + ray_test + ray:
                                    C_ent=False

                # on continue le parcours de la liste des inclusions existantes
                i = i + 1

            # On ajoute la nouvelle inclusion, si cela est possible ; on calcule aussi le nouveau volume occupe par les boules
            if C_ent:
                # ajout de l'inclusion, decalage entre le numero de phase et phi, choisi pur parcourir des tableaux
                liste_ph.append([cen, ray, phi+1])
                # calcul du nouveau volume occupe par la pĥase phi
                vol_inc[phi] = vol_inc[phi] + vol3D(ray)

        # fraction volumique occupee par les phases solides
        total_frac_vol = vol_inc/size**3

    return(liste_ph, total_frac_vol)
