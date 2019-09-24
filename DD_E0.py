###############################################################################################################
############################### Generation de structures periodiques aleatoires ###############################
###############################################################################################################

### ------------ Paquets a importer ------------ ###

# paquets mathematiques
import numpy as np
import random as rd
from math import sqrt
from math import exp

# affichage etc

import matplotlib

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gridspec

import sys, os
import shelve as sh

import multiprocessing


import pylab as pl
#from pylab import *

# cmap_conv = 'night'
cmap_conv = 'cement'

if cmap_conv == 'night':
    cmap = matplotlib.colors.ListedColormap(['blue','grey','yellow'])
elif cmap_conv == 'cement':
    cmap = matplotlib.colors.ListedColormap(['cyan','grey','orange'])

import time


# sys.exit('imports faits')

from RSAA_2d_ray import *

## Tache a accomplir : creer lise RSAA, matrice d'appartenance, representation graphique ##

# task = 'nothing'
task = 'listS'
# task = 'Ainc'
# task = 'graph'

## Generation d'inclusions dans un reseau carre, stockage ##

# cote=70
# Tps=10000

# cote=150
# Tps=20000

# cote=300
# Tps=25000

# cote=500
# Tps=30000

cote=750
Tps=45000

# cote=1000
# Tps=1000000

# geometrie euclidienne

geo=[eucl2D,vol2D]

# lois normales pour la taille des inclusions ; parametres #

Lois=[g_norm,g_norm]
Par=[(3,1),(1.5,5)]#[(5,0.5),(8,0.1)] pour 1000

# fractions volumiAincques pour les differentes phases ; marges

frac=[0.15,0.12]
delta_inc=2

n_proc = 40

rep_inc='IncAlea2D'









## --- Realisation des taches --- ##

start = time.time()

if task == 'listS':

    L = RSAA_ph_dist2D(geo,delta_inc,Lois,Par,frac,[cote,cote],Tps,'per')

    l_name='L_'+str(cote)+'_'+str(Tps)
    # stockage : listes de centres-rayons-phases et fractions volumiques #
    # sauvegarde de la variable maliste sous le nom "maliste" dans le fichier 'L_20_10000'
    with sh.open(rep_inc+'/'+l_name) as l_sto:
        l_sto["maliste"] = L

    print('liste RSAA faite, '+str(len(L[0]))+' inclusions')

elif task == 'Ainc':

    # chargement de la variable maliste, on recharge les variables de dimension spatiale et du nombre de tours de boucle RSAA

    #cote=
    #Tps=
    l_name = 'L_'+str(cote)+'_'+str(Tps)

    with sh.open(rep_inc+'/'+l_name) as l_loa:
        L_loa = l_loa["maliste"]

    # sys.exit('Remplissage de A, '+str(len(L_loa[0]))+' inclusions')
    print('Remplissage de A, '+str(len(L_loa[0]))+' inclusions')
    ## creation d'une matrice contenant les inclustions : un coefficient pour un pixel ##
    if n_proc == 1:
        A_remp = Vremp2D(L_loa,eucl2D,[cote,cote],'per')

    #pour avoir une fonction top-level e paralleliser

    L=L_loa[0]
    M=np.array(L)
    n_phi=max(M[:,2])

    dim=[cote,cote]
    C_per='per'
    dist=eucl2D

    def parc_liste_k(k):
        return(parc_liste(k,L,n_phi,dist,dim,C_per))

    # n_proc = 40

    # print('Processeurs',n_proc)

    pool=multiprocessing.Pool(processes=n_proc)

    B_par=pool.map(parc_liste_k,(k for k in range(dim[0]*dim[1])))

    A_remp=np.array(B_par)
    A_remp=np.reshape(A_remp,(dim[0],dim[1]))

    ##500 ; tps 30 000 ; 4600 inclusions : environ 1h sur MECALAC 8 x 2,90 GHz
    ##1000 ; tps 1 000 000 ; 2800 inclusions : environ 2 h 40 min pour l'ordinateur fixe 8 x 3,50GHz


    ##70 ; tps 10 000 ; ... inclusions : 20 secondes environ
    ##150 ; tps 20 000 ; 1000 inclusions : 2 minutes 20 secondes
    ##1000 ; tps 1 000 000 ; 2800 inclusions : 15h ++

    # stockage de la matrice A #

    a_name='VA_'+str(cote)+'_'+str(Tps)
    #repertoire ?

    rep_inc='IncAlea2D'

    with sh.open(rep_inc+'/'+a_name) as a_sto:
        a_sto["maliste"] = A_remp

    print('matrice A creee')

elif task == 'graph':

    a_name='VA_'+str(cote)+'_'+str(Tps)

    with sh.open(rep_inc+'/'+a_name) as a_loa:
        A_loa = a_loa["maliste"]

    ## representation graphique ##

    # Couleurs : bleu pour le fluide, gris pour une premiere inclusion, jaune pour une deuxieme
    ##avec pylab
    # cmap = matplotlib.colors.ListedColormap(['blue','grey','yellow'])
    bounds = [0,1,2,3]
    norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
    #

    pl.imshow(A_loa,interpolation='none',cmap=cmap,norm=norm)
    pl.axis('off')

    figname = cmap_conv + '_' + str(cote) + 'f' + str(cote) + 'Tps' + str(Tps) + '.png'
    rep = 'Figures2D'
    save_name = rep+'/'+figname

    pl.savefig(save_name)

    # pl.show()

end = time.time()

print('Tache effectuee, cote '+str(cote)+', '+str(end-start)+' secondes')
