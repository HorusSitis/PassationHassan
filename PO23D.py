# -*- coding: utf-8 -*-
"""
Created on Thu Jun  7 12:20:45 2018

@author: ghraieb
"""

# from fenics import *
# from math import *
# import numpy as np
# import matplotlib.pyplot as plt

#########################################
####Matrice de Correlation temporelle####
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
    ## on a effectue seulement num_steps*(num_step+1)/2 produits entre les fonctions ui-j

	for k in range(num_steps):
		ui.vector().set_local(U[:,k])
		S[k]=assemble(dot(ui,ui)*dx)
	S=S*(1.0/num_steps)
	N=np.diag(S)

	return C+N

##########################################################################################
#######################Matrice des coefficients aleatoires A (ai(t))######################
####Matrice de la base totale sans prendre en compte l'energie de chaque vecteur (phi)####
##########################################################################################

def mat_a_mat_phi(num_steps,U,C,V,base_POD_normee):
    ## Resolution du probleme aux valeurs-vecteurs propres, avec la matrice des correlations temporelles
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
    ## Resultats : valeurs propres ... ordre ? ; matrice de coefficients aleatoires A et base POD Phi_prime
    return [valp,A,Phi_prime]

##########################################################################################
##################### Energie et energie cumulee des valeurs propres #####################
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
##################### Calcul des coefficients du modele reduit, etant donne le domaine de definition du nouveau champ de vecteurs a calculer #####################
##################################################################################################################################################################

#V_nouv
#Phi_prime_v


tol=1e-10
def calc_Ab_2D(V_nouv,mesh_nouv,Phi_nouv_v,r_nouv,cen,nb_modes):
    A=np.zeros((nb_modes,nb_modes))
    b=np.zeros(nb_modes)
    ### Fonctions a definir pour calculer les coefficients des deux tenseurs, qui dependent de la metrique de l'espace des fonctions test
    phi_nouv_k=Function(V_nouv)
    phi_nouv_i=Function(V_nouv)
    # boucle pour le calcul de la matrice de coefficients
    for k in range(nb_modes):
        phi_nouv_k.vector().set_local(Phi_nouv_v[:,k])
        for i in range(nb_modes):
            phi_nouv_i.vector().set_local(Phi_nouv_v[:,i])
            # On calcule le coefficient Aki
            A[k,i]=assemble(tr(dot((grad(phi_nouv_k)).T, grad(phi_nouv_i)))*dx)
    # creation de l'interface solide-fluide
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
    # boucle pour le calcul du second membre du probleme lineaire MOR
    for i in range(nb_modes):
        phi_nouv_i.vector().set_local(Phi_nouv_v[:,i])
        b[i]=assemble(dot(normale,phi_nouv_i)*ds(num_front_inc))
    return([A,b])

def calc_Ab_3D(V_nouv,mesh_nouv,Phi_nouv_v,r_nouv,origin,nb_modes,config):
    A=np.zeros((nb_modes,nb_modes))
    b=np.zeros(nb_modes)
    ### Fonctions a definir pour calculer les coefficients des deux tenseurs, qui dependent de la metrique de l'espace des fonctions test
    phi_nouv_k=Function(V_nouv)
    phi_nouv_i=Function(V_nouv)
    # boucle pour le calcul de la matrice de coefficients
    for k in range(nb_modes):
        phi_nouv_k.vector().set_local(Phi_nouv_v[:,k])
        for i in range(nb_modes):
            phi_nouv_i.vector().set_local(Phi_nouv_v[:,i])
            # On calcule le coefficient Aki
            A[k,i]=assemble(tr(dot((grad(phi_nouv_k)).T, grad(phi_nouv_i)))*dx)
    # creation de l'interface solide-fluide
    r=r_nouv
    if config=='sph_un':
        l_cen=[origin]
        #print(l_cen)
        class inclusion_periodique(SubDomain):
            def inside(self,x,on_boundary):
                return (on_boundary and any([between((x[0]-c[0]), (-r-tol, r+tol)) for c in l_cen]) and any([between((x[1]-c[1]), (-r-tol, r+tol)) for c in l_cen]) and any([between((x[2]-c[2]), (-r-tol, r+tol)) for c in l_cen]))
    elif config=='cyl_un':
        l_axe=[origin]
        #print(l_axe)
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
    # boucle pour le calcul du second membre du probleme lineaire MOR
    for i in range(nb_modes):
        phi_nouv_i.vector().set_local(Phi_nouv_v[:,i])
        b[i]=assemble(dot(normale,phi_nouv_i)*ds(num_front_inc))
    return([A,b])

### Valable en dimension 2 ###

def calc_Ab_compl(V_nouv,mesh_nouv,Phi_nouv_v,nb_modes,test_snap):
    A=np.zeros((nb_modes,nb_modes))
    b=np.zeros(nb_modes)
    ## Fonctions a definir pour calculer les coefficients des deux tenseurs, qui dependent de la metrique de l'espace des fonctions test
    phi_nouv_k=Function(V_nouv)
    phi_nouv_i=Function(V_nouv)
    ## Boucle pour le calcul de la matrice de coefficients
    for k in range(nb_modes):
        phi_nouv_k.vector().set_local(Phi_nouv_v[:,k])
        for i in range(nb_modes):
            phi_nouv_i.vector().set_local(Phi_nouv_v[:,i])
            # On calcule le coefficient Aki
            A[k,i]=assemble(tr(dot((grad(phi_nouv_k)).T, grad(phi_nouv_i)))*dx)
    ## On definit la bordure du domaine, sur laquelle integrer le second membre "L" de l'equation en dimension finie
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
        print('Gamma sf ne coupe pas le bord du carre')
        boundaries.set_all(0)
        Gamma_sf.mark(boundaries, 1)
    elif test_snap=='solid_2':
        #
        num_front_inc=11
        print('Gamma sf coupe le bord du carre')
        #
    ## On integre les vecteurs POD pour obtenir les coefficients du modele reduit
    normale=FacetNormal(mesh_nouv)
    # boucle pour le calcul du second membre du probleme lineaire MOR
    for i in range(nb_modes):
        phi_nouv_i.vector().set_local(Phi_nouv_v[:,i])
        b[i]=assemble(dot(normale,phi_nouv_i)*ds(num_front_inc))
    return([A,b])

## Dimension 3 : meme algorithme ##

def calc_Ab_compl_3D(mesh_n_name,Phi_nouv_v,nb_modes):
    ## Maillage et espace de fonctions test, depuis le repertoire PassationHassan
    mesh_nouv=Mesh(mesh_n_name+".xml")
    V_nouv=VectorFunctionSpace(mesh_nouv, "P", 2, constrained_domain=PeriodicBoundary())
    ## Matrice de resultats, initialisees a 0
    A=np.zeros((nb_modes,nb_modes))
    b=np.zeros(nb_modes)
    ## Fonctions a definir pour calculer les coefficients des deux tenseurs, qui dependent de la metrique de l'espace des fonctions test
    phi_nouv_k=Function(V_nouv)
    phi_nouv_i=Function(V_nouv)
    ## Boucle pour le calcul de la matrice de coefficients
    for k in range(nb_modes):
        phi_nouv_k.vector().set_local(Phi_nouv_v[:,k])
        for i in range(nb_modes):
            phi_nouv_i.vector().set_local(Phi_nouv_v[:,i])
            # On calcule le coefficient Aki
            A[k,i]=assemble(tr(dot((grad(phi_nouv_k)).T, grad(phi_nouv_i)))*dx)
    ## On definit la bordure du domaine, sur laquelle integrer le second membre "L" de l'equation en dimension finie
    boundaries = MeshFunction('size_t', mesh_nouv, mesh_n_name+"_facet_region"+".xml")
    ds = Measure("ds")(subdomain_data=boundaries)
    ## Marquage des bordures pour la condition de Neumann
    num_front_inc=1700
    ## On integre les vecteurs POD pour obtenir les coefficients du modele reduit
    normale=FacetNormal(mesh_nouv)
    # boucle pour le calcul du second membre du probleme lineaire MOR
    for i in range(nb_modes):
        phi_nouv_i.vector().set_local(Phi_nouv_v[:,i])
        b[i]=assemble(dot(normale,phi_nouv_i)*ds(num_front_inc))
    return([A,b])

################################################################################################################################################################################################
## ------------------------------------------------------------------ Sans interpolation : directement sur le domaine fixe ------------------------------------------------------------------ ##
################################################################################################################################################################################################

def calc_Ab_simpl_2D_ninterpol(V_r_fixe,mesh_r_fixe,Phi_r_fixe_v,r_nouv,cen,nb_modes):
    A=np.zeros((nb_modes,nb_modes))
    b=np.zeros(nb_modes)
    r=r_nouv
    ### Integration des coefficients ROM sur le volume fluide : premier terme du probleme faible
    ## Sous-domaine de l'espace fixe correspondant au domaine fluide courant
    class DomPhysFluide(SubDomain):
        def inside(self, x, on_boundary):
            return True if ((x[0]-cen[0])**2+(x[1]-cen[1])**2>=r**2) else False
    # marquage du domaine
    dom_courant=DomPhysFluide()
    subdomains=MeshFunction('size_t',mesh_r_fixe,mesh_r_fixe.topology().dim())
    subdomains.set_all(1)
    dom_courant.mark(subdomains,12829)
    dxf=Measure("dx", domain=mesh_r_fixe, subdomain_data=subdomains)
    ## Fonctions a definir pour calculer les coefficients des deux tenseurs, qui dependent de la metrique de l'espace des fonctions test
    phi_r_fixe_k=Function(V_r_fixe)
    phi_r_fixe_i=Function(V_r_fixe)
    # boucle pour le calcul de la matrice de coefficients
    for k in range(nb_modes):
        phi_r_fixe_k.vector().set_local(Phi_r_fixe_v[:,k])
        for i in range(nb_modes):
            phi_r_fixe_i.vector().set_local(Phi_r_fixe_v[:,i])
            # On calcule le coefficient Aki
            A[k,i]=assemble(tr(dot((grad(phi_r_fixe_k)).T, grad(phi_r_fixe_i)))*dxf(12829))
    ### Integration des coefficients ROM sur l'interface solide-fluide : condition de Neumann pour le probleme faible
    ## Creation de l'interface solide-fluide
    class Front(SubDomain):
        def inside(self,x,on_boundary):
            return True if between((x[0]-cen[0])**2+(x[1]-cen[1])**2,(r**2+tol,r**2-tol)) else False
    Gamma_sf=Front()
    boundaries = MeshFunction("size_t", mesh_r_fixe, mesh_r_fixe.topology().dim()-1)
    boundaries.set_all(1)
    num_front_inc=7
    Gamma_sf.mark(boundaries,num_front_inc)
    ds = Measure("ds")(subdomain_data=boundaries)
    normale=FacetNormal(mesh_r_fixe)
    # boucle pour le calcul du second membre du probleme lineaire MOR
    for i in range(nb_modes):
        phi_r_fixe_i.vector().set_local(Phi_r_fixe_v[:,i])
        b[i]=assemble(dot(normale,phi_r_fixe_i)*ds(num_front_inc))
    c=assemble( dot( normale,normale )*ds(num_front_inc))
    return([A,b,c])




################################################################################################################################################################################################
################################################################################################################################################################################################





tol_sol=0.0001
def calc_Ab_simpl_3D_ninterpol(mesh_f_name,config,geo_p,r_nouv,Phi_prime_v,nb_modes):
    r=r_nouv
    ## Maillage et espace de fonctions test, depuis le repertoire PassationHassan
    mesh_fixe=Mesh(mesh_f_name)#+".xml")
    V_fixe=VectorFunctionSpace(mesh_fixe, "P", 2, constrained_domain=PeriodicBoundary())
    ## Domaine d'integration pour les coefficients de A
    if config=='sph_un' and geo_p=='ray':
        class DomPhysFluide(SubDomain):
            def inside(self, x, on_boundary):
                return True if ((x[0]-0.5)**2+(x[1]-0.5)**2+(x[2]-0.5)**2>=r**2) else False
        dom_courant=DomPhysFluide()
        subdomains=MeshFunction('size_t',mesh_fixe,mesh_fixe.topology().dim())
        subdomains.set_all(1)
        dom_courant.mark(subdomains,12829)
        dxf=Measure("dx", domain=mesh_fixe, subdomain_data=subdomains)
    elif config=='cyl_un' and geo_p=='ray':
        class DomPhysFluide(SubDomain):
            def inside(self, x, on_boundary):
                return True if ((x[0]-0.5)**2+(x[2]-0.5)**2>=r**2) else False
        dom_courant=DomPhysFluide()
        subdomains=MeshFunction('size_t',mesh_fixe,mesh_fixe.topology().dim())
        subdomains.set_all(1)
        dom_courant.mark(subdomains,12829)
        dxf=Measure("dx", domain=mesh_fixe, subdomain_data=subdomains)
    ## Creation de l'interface solide-fluide
    if config=='sph_un':
        cen=[0.5,0.5,0.5]
        class frontiere_physique(SubDomain):
            def inside(self,x,on_boundary):
                return between((x[0]-cen[0])**2+(x[1]-cen[1])**2+(x[2]-cen[2])**2,(r**2-tol_sol,r**2+tol_sol))
    elif config=='cyl_un':
        top=[0.5,0.,0.5]
        #print(l_axe)
        class frontiere_physique(SubDomain):
            def inside(self,x,on_boundary):
                return between((x[0]-top[0])**2+(x[2]-top[2])**2,(r**2-tol_sol,r**2+tol_sol))
    Gamma_sf=frontiere_physique()
    boundaries = MeshFunction("size_t", mesh_fixe, mesh_fixe.topology().dim()-1)
    boundaries.set_all(1)
    Gamma_sf.mark(boundaries,1)# 7)
    ds = Measure("ds")(subdomain_data=boundaries)
    num_ff=1
    num_front_inc=1#7
    normale=FacetNormal(mesh_fixe)
    ## Matrice de resultats, initialisees a 0
    A=np.zeros((nb_modes,nb_modes))
    b=np.zeros(nb_modes)
    ## Fonctions a definir pour calculer les coefficients des deux tenseurs, qui dependent de la metrique de l'espace des fonctions test
    phi_prime_k=Function(V_fixe)
    phi_prime_i=Function(V_fixe)
    ## Boucle pour le calcul de la matrice de coefficients
    for k in range(nb_modes):
        phi_prime_k.vector().set_local(Phi_prime_v[:,k])
        for i in range(nb_modes):
            phi_prime_i.vector().set_local(Phi_prime_v[:,i])
            # On calcule le coefficient Aki
            A[k,i]=assemble(tr(dot((grad(phi_prime_k)).T, grad(phi_prime_i)))*dxf(12829))
    ## Boucle pour le calcul du second membre du probleme lineaire MOR
    for i in range(nb_modes):
        phi_prime_i.vector().set_local(Phi_prime_v[:,i])
        b[i]=assemble(dot(normale,phi_prime_i)*ds(num_front_inc))
    return([A,b])





def calc_Ab_compl_3D_ninterpol(mesh_f_name,config,geo_p,r_cen,r_per,Phi_prime_v,nb_modes):
    ## Maillage et espace de fonctions test, depuis le repertoire PassationHassan
    mesh_fixe=Mesh(mesh_f_name)#+".xml")
    V_fixe=VectorFunctionSpace(mesh_fixe, "P", 2, constrained_domain=PeriodicBoundary())
    ## Domaine d'integration pour les coefficients de A
    if config=='2sph':
        class DomPhysFluide(SubDomain):
            def inside(self, x, on_boundary):
                return True if ((x[0]-0.5)**2+(x[1]-0.5)**2+(x[2]-0.5)**2>=r_cen**2) else False
        dom_courant=DomPhysFluide()
        subdomains=MeshFunction('size_t',mesh_fixe,mesh_fixe.topology().dim())
        subdomains.set_all(1)
        dom_courant.mark(subdomains,12829)
        dxf=Measure("dx", domain=mesh_fixe, subdomain_data=subdomains)
    elif config=='cylsph':
        class DomPhysFluide(SubDomain):
            def inside(self, x, on_boundary):
                return True if ((geo_p=='ray_sph' and ((x[0]-0.5)**2+(x[1]-0.5)**2+(x[2]-0.5)**2>=r_cen**2)) or (geo_p=='ray_cyl' and x[0]**2+x[2]**2>=r_per**2 and (1-x[0])**2+x[2]**2>=r_per**2 and x[0]**2+(1-x[2])**2>=r_per**2 and (1-x[0])**2+(1-x[2])**2>=r_per**2)) else False
        dom_courant=DomPhysFluide()
        subdomains=MeshFunction('size_t',mesh_fixe,mesh_fixe.topology().dim())
        subdomains.set_all(1)
        dom_courant.mark(subdomains,12829)
        dxf=Measure("dx", domain=mesh_fixe, subdomain_data=subdomains)
    ## Domaine d'integration pour les coefficients de b : condition de Neumann
    boundaries = MeshFunction('size_t', mesh_fixe, mesh_n_name+"_facet_region"+".xml")
    ds = Measure("ds")(subdomain_data=boundaries)
    # Marquage des bordures pour la condition de Neumann
    num_front_inc=1700
    # On integre les vecteurs POD pour obtenir les coefficients du modele reduit
    normale=FacetNormal(mesh_fixe)
    ## Matrice de resultats, initialisees a 0
    A=np.zeros((nb_modes,nb_modes))
    b=np.zeros(nb_modes)
    ## Fonctions a definir pour calculer les coefficients des deux tenseurs, qui dependent de la metrique de l'espace des fonctions test
    phi_prime_k=Function(V_fixe)
    phi_prime_i=Function(V_fixe)
    ## Boucle pour le calcul de la matrice de coefficients
    for k in range(nb_modes):
        phi_prime_k.vector().set_local(Phi_prime_v[:,k])
        for i in range(nb_modes):
            phi_prime_i.vector().set_local(Phi_prime_v[:,i])
            # On calcule le coefficient Aki
            A[k,i]=assemble(tr(dot((grad(phi_prime_k)).T, grad(phi_prime_i)))*dxf(12829))
    ## Boucle pour le calcul du second membre du probleme lineaire MOR
    for i in range(nb_modes):
        phi_prime_i.vector().set_local(Phi_prime_v[:,i])
        b[i]=assemble(dot(normale,phi_prime_i)*ds(num_front_inc))
    return([A,b])
