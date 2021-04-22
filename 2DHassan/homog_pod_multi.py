# -*- coding: utf-8 -*-
"""
Created on Wed Sep  5 10:44:23 2018

@author: ghraieb
"""

from fenics import *
from dolfin import *
### --- from mshr import * --- ###
import matplotlib.pyplot as plt
import numpy as np
from calc_POD import *
from Lecture_ecriture_homog import *
from math import sqrt
import sys 

tol=1e-10

xinf=0.0
yinf=0.0
xsup=1.0
ysup=1.0

c_x=0.5
c_y=0.5

c1_x=0.0
c1_y=0.0

c2_x=1.0
c2_y=0.0

c3_x=1.0
c3_y=1.0

c4_x=0.0
c4_y=1.0

R_1=np.arange(0.1,0.5,0.05)
R_2=np.arange(0.05,0.5,0.05)
R_3=[0.01,0.05,0.1,0.15,0.2,0.22,0.24,0.25,0.255]

R_dim1=len(R_1)
R_dim2=len(R_2)
R_dim3=len(R_3)
R_dim=R_dim1+R_dim2+R_dim3


def raffinemment_maillage_1(r,mesh):
    markers = CellFunction("bool", mesh)
    markers.set_all(False)
    for c in cells(mesh):
	# Mark cells with facet midpoints near y == 1.0
	for f in facets(c):
	    if (sqrt((f.midpoint()[0]-c_x)**2+(f.midpoint()[1]-c_y)**2)<=1.2*r):
		markers[c] = True
    #new_mesh=refine(mesh, markers, redistribute=False)  
    new_mesh=refine(mesh, markers, redistribute=True)  
    
    return new_mesh

def raffinemment_maillage_2(r,mesh):
    markers = CellFunction("bool", mesh)
    markers.set_all(False)
    for c in cells(mesh):
	# Mark cells with facet midpoints near y == 1.0
	for f in facets(c):
	    if ((sqrt((f.midpoint()[0]-c1_x)**2+(f.midpoint()[1]-c1_y)**2)<=1.2*r) or (sqrt((f.midpoint()[0]-c2_x)**2+(f.midpoint()[1]-c2_y)**2)<=1.2*r) or (sqrt((f.midpoint()[0]-c3_x)**2+(f.midpoint()[1]-c3_y)**2)<=1.2*r) or (sqrt((f.midpoint()[0]-c4_x)**2+(f.midpoint()[1]-c4_y)**2)<=1.2*r)):
		markers[c] = True
    #new_mesh=refine(mesh, markers, redistribute=False)  
    new_mesh=refine(mesh, markers, redistribute=True)  
    
    return new_mesh

def raffinemment_maillage_3(r,mesh):
    markers = CellFunction("bool", mesh)
    markers.set_all(False)
    for c in cells(mesh):
	# Mark cells with facet midpoints near y == 1.0
	for f in facets(c):
	    if ((sqrt((f.midpoint()[0]-c_x)**2+(f.midpoint()[1]-c_y)**2)<=1.2*0.45) or (sqrt((f.midpoint()[0]-c1_x)**2+(f.midpoint()[1]-c1_y)**2)<=1.2*r) or (sqrt((f.midpoint()[0]-c2_x)**2+(f.midpoint()[1]-c2_y)**2)<=1.2*r) or (sqrt((f.midpoint()[0]-c3_x)**2+(f.midpoint()[1]-c3_y)**2)<=1.2*r) or (sqrt((f.midpoint()[0]-c4_x)**2+(f.midpoint()[1]-c4_y)**2)<=1.2*r)):
		markers[c] = True
    #new_mesh=refine(mesh, markers, redistribute=False)  
    new_mesh=refine(mesh, markers, redistribute=True)  
    
    return new_mesh

def raffinemment_maillage_q_1(r,mesh):
    markers = CellFunction("bool", mesh)
    markers.set_all(False)
    for c in cells(mesh):
	# Mark cells with facet midpoints near y == 1.0
	for f in facets(c):
	    if ((sqrt((f.midpoint()[0]-c_x)**2+(f.midpoint()[1]-c_y)**2)<=1.2*0.45) or (sqrt((f.midpoint()[0]-c1_x)**2+(f.midpoint()[1]-c1_y)**2)<=1.2*r)):
		markers[c] = True
    #new_mesh=refine(mesh, markers, redistribute=False)  
    new_mesh=refine(mesh, markers, redistribute=True)  
    
    return new_mesh

def raffinemment_maillage_q_2(r,mesh):
    markers = CellFunction("bool", mesh)
    markers.set_all(False)
    for c in cells(mesh):
	# Mark cells with facet midpoints near y == 1.0
	for f in facets(c):
	    if ((sqrt((f.midpoint()[0]-c_x)**2+(f.midpoint()[1]-c_y)**2)<=1.2*0.45) or (sqrt((f.midpoint()[0]-c1_x)**2+(f.midpoint()[1]-c1_y)**2)<=1.2*r) or (sqrt((f.midpoint()[0]-c2_x)**2+(f.midpoint()[1]-c2_y)**2)<=1.2*r) or (sqrt((f.midpoint()[0]-c3_x)**2+(f.midpoint()[1]-c3_y)**2)<=1.2*r)):
		markers[c] = True
    #new_mesh=refine(mesh, markers, redistribute=False)  
    new_mesh=refine(mesh, markers, redistribute=True)  
    
    return new_mesh

def raffinemment_maillage_q_3(r,mesh):
    markers = CellFunction("bool", mesh)
    markers.set_all(False)
    for c in cells(mesh):
	# Mark cells with facet midpoints near y == 1.0
	for f in facets(c):
	    if ((sqrt((f.midpoint()[0]-c_x)**2+(f.midpoint()[1]-c_y)**2)<=1.2*0.45) or (sqrt((f.midpoint()[0]-c1_x)**2+(f.midpoint()[1]-c1_y)**2)<=1.2*r) or (sqrt((f.midpoint()[0]-c2_x)**2+(f.midpoint()[1]-c2_y)**2)<=1.2*r)):
		markers[c] = True
    #new_mesh=refine(mesh, markers, redistribute=False)  
    new_mesh=refine(mesh, markers, redistribute=True)  
    
    return new_mesh

def raffinemment_maillage_q_4(r,mesh):
    markers = CellFunction("bool", mesh)
    markers.set_all(False)
    for c in cells(mesh):
	# Mark cells with facet midpoints near y == 1.0
	for f in facets(c):
	    if ((sqrt((f.midpoint()[0]-c_x)**2+(f.midpoint()[1]-c_y)**2)<=1.2*0.45) or (sqrt((f.midpoint()[0]-c1_x)**2+(f.midpoint()[1]-c1_y)**2)<=1.2*r) or (sqrt((f.midpoint()[0]-c3_x)**2+(f.midpoint()[1]-c3_y)**2)<=1.2*r)):
		markers[c] = True
    #new_mesh=refine(mesh, markers, redistribute=False)  
    new_mesh=refine(mesh, markers, redistribute=True)  
    
    return new_mesh

def raffinemment_maillage_q_5(r1,r2,r3,r4,mesh):
    markers = CellFunction("bool", mesh)
    markers.set_all(False)
    for c in cells(mesh):
	# Mark cells with facet midpoints near y == 1.0
	for f in facets(c):
	    if ((sqrt((f.midpoint()[0]-c_x)**2+(f.midpoint()[1]-c_y)**2)<=1.2*0.45) or (sqrt((f.midpoint()[0]-c1_x)**2+(f.midpoint()[1]-c1_y)**2)<=1.2*r1) or (sqrt((f.midpoint()[0]-c2_x)**2+(f.midpoint()[1]-c2_y)**2)<=1.2*r2) or (sqrt((f.midpoint()[0]-c3_x)**2+(f.midpoint()[1]-c3_y)**2)<=1.2*r3) or (sqrt((f.midpoint()[0]-c4_x)**2+(f.midpoint()[1]-c4_y)**2)<=1.2*r4)):
		markers[c] = True
    #new_mesh=refine(mesh, markers, redistribute=False)  
    new_mesh=refine(mesh, markers, redistribute=True)  
    
    return new_mesh

def raffinemment_maillage_q_6(r1,r2,r3,r4,mesh):
    markers = CellFunction("bool", mesh)
    markers.set_all(False)
    for c in cells(mesh):
	# Mark cells with facet midpoints near y == 1.0
	for f in facets(c):
	    if ((sqrt((f.midpoint()[0]-c1_x)**2+(f.midpoint()[1]-c1_y)**2)<=1.2*r1) or (sqrt((f.midpoint()[0]-c2_x)**2+(f.midpoint()[1]-c2_y)**2)<=1.2*r2) or (sqrt((f.midpoint()[0]-c3_x)**2+(f.midpoint()[1]-c3_y)**2)<=1.2*r3) or (sqrt((f.midpoint()[0]-c4_x)**2+(f.midpoint()[1]-c4_y)**2)<=1.2*r4)):
		markers[c] = True
    #new_mesh=refine(mesh, markers, redistribute=False)  
    new_mesh=refine(mesh, markers, redistribute=True)  
    
    return new_mesh


# Sub domain for Periodic boundary condition
class PeriodicBoundary(SubDomain):
        
        # Left boundary is "target domain" G
    def inside(self, x, on_boundary):
        return on_boundary and (near(x[0],xinf,tol) or near(x[1],yinf,tol))

            
            # Map right boundary (H) to left boundary (G)
    def map(self, x, y):
        if (near(x[0],xsup,tol)):
            y[0] = x[0] - 1.0
            y[1] = x[1]
                
        else :
            y[0]=x[0]
            y[1] = x[1] - 1.0


domaine_fixe=Rectangle(Point(xinf,yinf),Point(xsup,ysup))
mesh_fixe=generate_mesh(domaine_fixe,100)

V_fixe=VectorFunctionSpace(mesh_fixe, "P", 2, constrained_domain=PeriodicBoundary())
nb_noeuds = V_fixe.dim()


Usnap=np.zeros((nb_noeuds,R_dim))

repertoire_read_1="Solutions_homog_interp_circulaire/"
repertoire_read_2="Solutions_homog_interp_complexe_1/"
repertoire_read_3="Solutions_homog_interp_complexe/"
#repertoire_read_q="Solutions_homog_interp_quelconque/"


u = Function(V_fixe)

for n in range(R_dim1):
    r=R_1[n]
    lecture_fichiers_restart_hdf5(u,repertoire_read_1,r,mesh_fixe) 
    Usnap[:,n]=u.vector().array()

for n in range(R_dim2):
    r=R_2[n]
    lecture_fichiers_restart_hdf5(u,repertoire_read_2,r,mesh_fixe) 
    Usnap[:,R_dim1+n]=u.vector().array()

for n in range(R_dim3):
    r=R_3[n]
    lecture_fichiers_restart_hdf5(u,repertoire_read_3,r,mesh_fixe) 
    Usnap[:,R_dim1+R_dim2+n]=u.vector().array()

#matrice de corrélation
C=mat_corr_temp(V_fixe,R_dim,Usnap)

# Calcul des coefficients aléatoires et la base POD
val_propres,A,phi=mat_a_mat_phi(R_dim,Usnap,C)

# Choix du nb de modes #
####pourcentage d'energie par val_propres#####
nb_modes=0
ener_pour=np.zeros((R_dim))
ener_pour_cumul=np.zeros((R_dim))
s_t=0
for k in range(R_dim):
    s_t=s_t+val_propres[k]
for i in range(R_dim):
    s=val_propres[i]
    s_n=0
    for j in range(i+1):
        s_n=s_n+val_propres[j]
    ener_pour[i]=(s/s_t)*100
    ener_pour_cumul[i]=(s_n/s_t)*100
for i in range(R_dim):
    if (ener_pour_cumul[i] >= 99.99):
        nb_modes=i+1
        break
absc=np.arange(1,R_dim+1,1)
plt.plot(absc,ener_pour)
plt.show()
plt.plot(absc,ener_pour_cumul)
plt.show()

nb_modes=25

print nb_modes

#choix de la base POD
G=np.zeros((nb_noeuds,nb_modes))
for i in range(nb_modes):
    G[:,i]=phi[:,i]




#################################################################
###########################circulaire############################
#################################################################
#########Interpolation de la base dans chaque maillage###########
#################################################################
#
#
#r=0.423
#
#rect= Rectangle(Point(xinf,yinf),Point(xsup,ysup))
#circle = Circle(Point(c_x,c_y),r)
#domain=rect-circle
#res = 50  # Resolution of mesh
#mesh = generate_mesh(domain, res)
#
##### Premier raffinement dans l axe horizontale du cylindre 
#meshB= raffinemment_maillage_1(r,mesh)
#mesh=meshB
##plot(mesh)
##plt.show()
#   
#
##class Obstacle(SubDomain):
##    def inside(self, x, on_boundary):
##        return (on_boundary and between((x[0]-c_x), (-r-tol, r+tol)) and between((x[1]-c_x), (-r-tol, r+tol)))
#class Obstacle(SubDomain):
#    def inside(self, x, on_boundary):
#        return (on_boundary and between((x[0]-c_x), (-r-tol, r+tol)) and between((x[1]-c_x), (-r-tol, r+tol)))
#
#
#V=VectorFunctionSpace(mesh, 'P', 2, constrained_domain=PeriodicBoundary())
#nb_noeuds_ess=V.dim()
#
#print ('uuuuuuu',nb_noeuds_ess)
#
#obstacle = Obstacle()
#    
#P_impose=Constant((0., 0.))
#pref = DirichletBC(V, P_impose, "x[0] < DOLFIN_EPS && x[1] < DOLFIN_EPS", "pointwise")
#    
#
#
#boundaries = FacetFunction("size_t", mesh)
#
#boundaries.set_all(0)
#obstacle.mark(boundaries, 5)
#
#ds = Measure("ds")[boundaries]
#normale = FacetNormal(mesh) ## normale associe a chaque cellule
#    
#num_front_cercle=5 ## numero de la frontiere du cylindre
#
##################################################
######Sauvegarder les phi apres interpolation#####
##################################################  
#w=np.zeros((nb_noeuds_ess,nb_modes))
#wi_fixe=Function(V_fixe)
#
#repertoire_final_1="Homog_base_maill_init_multi/"
#
#File("Homog_base_maill_init_multi/mesh_circulaire.xml.gz") << mesh
#
#
#ufile,USAVE=creation_fichier_pourecriture_champ_hdf5(repertoire_final_1,mesh)
#file_rayon_ecriture = open("%s/rayon_ecriture.txt" %(repertoire_final_1), "w")
#kfic=0
#for n in range(nb_modes):
#    xc=n+1
#    wi_fixe.vector().set_local(G[:,n])
#    wi_fixe.set_allow_extrapolation(True)
#    wi_ess=interpolate(wi_fixe,V)
##    plot(wi_ess)
##    plt.show()
#    ecriture_champ_hdf5(ufile,USAVE,wi_ess,kfic,file_rayon_ecriture,xc)
#    kfic=kfic+1
#    w[:,n]=wi_ess.vector().array()
#    
#
#
######Model reduit pour un rayon précis#####
#mat_mass=np.zeros((nb_modes,nb_modes))
#vec_sf=np.zeros((nb_modes))
#phiu_i = Function(V)                                                                    
#phiu_j = Function(V)                                                                    
#    
#
#
#for iproj in range(nb_modes):
#    phiu_i.vector().set_local(w[:,iproj])
#	
#        
#    for jrec in range(nb_modes):
#          phiu_j.vector().set_local(w[:,jrec])
#              
#              
#          ### matrice de masse
#          mat_mass[jrec,iproj]=assemble(tr(dot((grad(phiu_j)).T, grad(phiu_i)))*dx)
#          
#
#for iproj in range(nb_modes):
#    phiu_i.vector().set_local(w[:,iproj])
#    ### terme sf
#    vec_sf[iproj]=assemble(dot(normale,phiu_i)*ds(num_front_cercle))
#
#sol=np.linalg.solve(mat_mass.T,-vec_sf)
#
#print sol
#
#repertoire_final_2="Solution_homog_red_maill_init_multi/"
#
#ufile,USAVE=creation_fichier_pourecriture_champ_hdf5(repertoire_final_2,mesh)
#file_rayon_ecriture = open("%s/rayon_ecriture.txt" %(repertoire_final_2), "w")
#kfic=0
#
#
#urec=np.dot(w,sol)
#urec_red=Function(V)
#urec_red.vector().set_local(urec)
#
#ecriture_champ_hdf5(ufile,USAVE,urec_red,kfic,file_rayon_ecriture,1)
#
#H=assemble(grad(urec_red)[0,0]*dx)
#M=assemble(grad(urec_red)[1,1]*dx)
#X=assemble(grad(urec_red)[0,1]*dx)
#C=assemble(grad(urec_red)[1,0]*dx)
#print H,M,X,C
#ep=(1-pi*(r**2))
#D=((H/ep)+1)*ep
#print ("hi",r,ep,D)
#
#
#
#plot(urec_red)
#plt.show()



#################################################################
##########################complexe_1#############################
#################################################################
#########Interpolation de la base dans chaque maillage###########
#################################################################
#
#
#r=0.087
#
#rect= Rectangle(Point(xinf,yinf),Point(xsup,ysup))
#circle1 = Circle(Point(c1_x,c1_y),r)
#circle2 = Circle(Point(c2_x,c2_y),r)
#circle3 = Circle(Point(c3_x,c3_y),r)
#circle4 = Circle(Point(c4_x,c4_y),r)
#domain=rect-(circle1+circle2+circle3+circle4)
#res = 50  # Resolution of mesh
#mesh = generate_mesh(domain, res)
#
##### Premier raffinement dans l axe horizontale du cyٍlindre 
#meshB= raffinemment_maillage_2(r,mesh)
#mesh=meshB
#
#
#class Obstacle1(SubDomain):
#    def inside(self, x, on_boundary):
#        return (on_boundary and between((x[0]-c1_x), (-r-tol, r+tol)) and between((x[1]-c1_y), (-r-tol, r+tol)))
#
#
#class Obstacle2(SubDomain):
#    def inside(self, x, on_boundary):
#        return (on_boundary and between((x[0]-c2_x), (-r-tol, r+tol)) and between((x[1]-c2_y), (-r-tol, r+tol)))
#
#
#class Obstacle3(SubDomain):
#    def inside(self, x, on_boundary):
#        return (on_boundary and between((x[0]-c3_x), (-r-tol, r+tol)) and between((x[1]-c3_y), (-r-tol, r+tol)))
#
#
#class Obstacle4(SubDomain):
#    def inside(self, x, on_boundary):
#        return (on_boundary and between((x[0]-c4_x), (-r-tol, r+tol)) and between((x[1]-c4_y), (-r-tol, r+tol)))
#
#
#
#V=VectorFunctionSpace(mesh, 'P', 2, constrained_domain=PeriodicBoundary())
#nb_noeuds_ess=V.dim()
#
#print ('uuuuuuu',nb_noeuds_ess)
#
#obstacle1 = Obstacle1()
#obstacle2 = Obstacle2()
#obstacle3 = Obstacle3()
#obstacle4 = Obstacle4()
#
#boundaries = FacetFunction("size_t", mesh)
#
#boundaries.set_all(0)
#
#obstacle1.mark(boundaries, 1)
#obstacle2.mark(boundaries, 2)
#obstacle3.mark(boundaries, 3)
#obstacle4.mark(boundaries, 4)
#
#    
#ds = Measure("ds")[boundaries]
#normale = FacetNormal(mesh) ## normale associe a chaque cellule
#
#
#num1=1
#num2=2
#num3=3
#num4=4
#
#
##################################################
######Sauvegarder les phi apres interpolation#####
##################################################  
#w=np.zeros((nb_noeuds_ess,nb_modes))
#wi_fixe=Function(V_fixe)
#
#repertoire_final_1="Homog_base_maill_init_multi/"
#
#File("Homog_base_maill_init_multi/mesh_complexe_1.xml.gz") << mesh
#
#ufile,USAVE=creation_fichier_pourecriture_champ_hdf5(repertoire_final_1,mesh)
#file_rayon_ecriture = open("%s/rayon_ecriture.txt" %(repertoire_final_1), "w")
#kfic=0
#for n in range(nb_modes):
#    xc=n+1
#    wi_fixe.vector().set_local(G[:,n])
#    wi_fixe.set_allow_extrapolation(True)
#    wi_ess=interpolate(wi_fixe,V)
##    plot(wi_ess)
##    plt.show()
#    ecriture_champ_hdf5(ufile,USAVE,wi_ess,kfic,file_rayon_ecriture,xc)
#    kfic=kfic+1
#    w[:,n]=wi_ess.vector().array()
#    
#
#
######Model reduit pour un rayon précis#####
#mat_mass=np.zeros((nb_modes,nb_modes))
#vec_sf=np.zeros((nb_modes))
#phiu_i = Function(V)                                                                    
#phiu_j = Function(V)                                                                    
#    
#
#
#for iproj in range(nb_modes):
#    phiu_i.vector().set_local(w[:,iproj])
#	
#        
#    for jrec in range(nb_modes):
#          phiu_j.vector().set_local(w[:,jrec])
#              
#              
#          ### matrice de masse
#          mat_mass[jrec,iproj]=assemble(tr(dot((grad(phiu_j)).T, grad(phiu_i)))*dx)
#          
#
#for iproj in range(nb_modes):
#    phiu_i.vector().set_local(w[:,iproj])
#    ### terme sf
#    vec_sf[iproj]=assemble((dot(normale,phiu_i)*ds(num1))+(dot(normale,phiu_i)*ds(num2))+(dot(normale,phiu_i)*ds(num3))+(dot(normale,phiu_i)*ds(num4)))
#
#sol=np.linalg.solve(mat_mass.T,-vec_sf)
#
#print sol
#
#repertoire_final_2="Solution_homog_red_maill_init_multi/"
#
#ufile,USAVE=creation_fichier_pourecriture_champ_hdf5(repertoire_final_2,mesh)
#file_rayon_ecriture = open("%s/rayon_ecriture.txt" %(repertoire_final_2), "w")
#kfic=0
#
#
#urec=np.dot(w,sol)
#urec_red=Function(V)
#urec_red.vector().set_local(urec)
#
#ecriture_champ_hdf5(ufile,USAVE,urec_red,kfic,file_rayon_ecriture,1)
#
#H=assemble(grad(urec_red)[0,0]*dx)
#M=assemble(grad(urec_red)[1,1]*dx)
#X=assemble(grad(urec_red)[0,1]*dx)
#C=assemble(grad(urec_red)[1,0]*dx)
#print H,M,X,C
#ep=(1-pi*(r**2))
#D=((H/ep)+1)*ep
#print ("hi",r,ep,D)
#
##print urec_red(0,0.5),urec_red(1,0.5),urec_red(0.5,0),urec_red(0.5,1),urec_red(0,0.4),urec_red(1,0.4)
#
#plot(urec_red)
#plt.show()




#################################################################
##############################complexe###########################
#################################################################
#########Interpolation de la base dans chaque maillage###########
#################################################################
#
#
#r=0.248
#
#rect= Rectangle(Point(xinf,yinf),Point(xsup,ysup))
#circle = Circle(Point(c_x,c_y),0.45)
#circle1 = Circle(Point(c1_x,c1_y),r)
#circle2 = Circle(Point(c2_x,c2_y),r)
#circle3 = Circle(Point(c3_x,c3_y),r)
#circle4 = Circle(Point(c4_x,c4_y),r)
#domain=rect-(circle+circle1+circle2+circle3+circle4)
#res = 50  # Resolution of mesh
#mesh = generate_mesh(domain, res)
#
##### Premier raffinement dans l axe horizontale du cyٍlindre 
#meshB= raffinemment_maillage_3(r,mesh)
#mesh=meshB
#
#class Obstacle(SubDomain):
#    def inside(self, x, on_boundary):
#        return (on_boundary and between((x[0]-c_x), (-0.45-tol, 0.45+tol)) and between((x[1]-c_x), (-0.45-tol, 0.45+tol)))
#
#
#class Obstacle1(SubDomain):
#    def inside(self, x, on_boundary):
#        return (on_boundary and between((x[0]-c1_x), (-r-tol, r+tol)) and between((x[1]-c1_y), (-r-tol, r+tol)))
#
#
#class Obstacle2(SubDomain):
#    def inside(self, x, on_boundary):
#        return (on_boundary and between((x[0]-c2_x), (-r-tol, r+tol)) and between((x[1]-c2_y), (-r-tol, r+tol)))
#
#
#class Obstacle3(SubDomain):
#    def inside(self, x, on_boundary):
#        return (on_boundary and between((x[0]-c3_x), (-r-tol, r+tol)) and between((x[1]-c3_y), (-r-tol, r+tol)))
#
#
#class Obstacle4(SubDomain):
#    def inside(self, x, on_boundary):
#        return (on_boundary and between((x[0]-c4_x), (-r-tol, r+tol)) and between((x[1]-c4_y), (-r-tol, r+tol)))
#
#
#V=VectorFunctionSpace(mesh, 'P', 2, constrained_domain=PeriodicBoundary())
#nb_noeuds_ess=V.dim()
#
#print ('uuuuuuu',nb_noeuds_ess)
#
#obstacle = Obstacle()
#obstacle1 = Obstacle1()
#obstacle2 = Obstacle2()
#obstacle3 = Obstacle3()
#obstacle4 = Obstacle4()
#
#boundaries = FacetFunction("size_t", mesh)
#
#boundaries.set_all(5)
#obstacle.mark(boundaries, 0)
#obstacle1.mark(boundaries, 1)
#obstacle2.mark(boundaries, 2)
#obstacle3.mark(boundaries, 3)
#obstacle4.mark(boundaries, 4)
#
#    
#ds = Measure("ds")[boundaries]
#normale = FacetNormal(mesh) ## normale associe a chaque cellule
#
#num0=0
#num1=1
#num2=2
#num3=3
#num4=4
#
#
##################################################
######Sauvegarder les phi apres interpolation#####
##################################################  
#w=np.zeros((nb_noeuds_ess,nb_modes))
#wi_fixe=Function(V_fixe)
#
#repertoire_final_1="Homog_base_maill_init_multi/"
#
#File("Homog_base_maill_init_multi/mesh_complexe.xml.gz") << mesh
#
#ufile,USAVE=creation_fichier_pourecriture_champ_hdf5(repertoire_final_1,mesh)
#file_rayon_ecriture = open("%s/rayon_ecriture.txt" %(repertoire_final_1), "w")
#kfic=0
#for n in range(nb_modes):
#    xc=n+1
#    wi_fixe.vector().set_local(G[:,n])
#    wi_fixe.set_allow_extrapolation(True)
#    wi_ess=interpolate(wi_fixe,V)
##    plot(wi_ess)
##    plt.show()
#    ecriture_champ_hdf5(ufile,USAVE,wi_ess,kfic,file_rayon_ecriture,xc)
#    kfic=kfic+1
#    w[:,n]=wi_ess.vector().array()
#    
#
#
######Model reduit pour un rayon précis#####
#mat_mass=np.zeros((nb_modes,nb_modes))
#vec_sf=np.zeros((nb_modes))
#phiu_i = Function(V)                                                                    
#phiu_j = Function(V)                                                                    
#    
#
#
#for iproj in range(nb_modes):
#    phiu_i.vector().set_local(w[:,iproj])
#	
#        
#    for jrec in range(nb_modes):
#          phiu_j.vector().set_local(w[:,jrec])
#              
#              
#          ### matrice de masse
#          mat_mass[jrec,iproj]=assemble(tr(dot((grad(phiu_j)).T, grad(phiu_i)))*dx)
#          
#
#for iproj in range(nb_modes):
#    phiu_i.vector().set_local(w[:,iproj])
#    ### terme sf
#    vec_sf[iproj]=assemble((dot(normale,phiu_i)*ds(num0))+(dot(normale,phiu_i)*ds(num1))+(dot(normale,phiu_i)*ds(num2))+(dot(normale,phiu_i)*ds(num3))+(dot(normale,phiu_i)*ds(num4)))
#
#sol=np.linalg.solve(mat_mass.T,-vec_sf)
#
#print sol
#
#repertoire_final_2="Solution_homog_red_maill_init_multi/"
#
#ufile,USAVE=creation_fichier_pourecriture_champ_hdf5(repertoire_final_2,mesh)
#file_rayon_ecriture = open("%s/rayon_ecriture.txt" %(repertoire_final_2), "w")
#kfic=0
#
#
#urec=np.dot(w,sol)
#urec_red=Function(V)
#urec_red.vector().set_local(urec)
#
#ecriture_champ_hdf5(ufile,USAVE,urec_red,kfic,file_rayon_ecriture,1)
#
#H=assemble(grad(urec_red)[0,0]*dx)
#M=assemble(grad(urec_red)[1,1]*dx)
#X=assemble(grad(urec_red)[0,1]*dx)
#C=assemble(grad(urec_red)[1,0]*dx)
#print H,M,X,C
#ep=(1-pi*((r**2)+0.45**2))
#D=((H/ep)+1)*ep
#print ("hi",r,ep,D)
#
#
#
#plot(urec_red)
#plt.show()





#################################################################
###################un quart et cercle central####################
#################################################################
#########Interpolation de la base dans chaque maillage###########
#################################################################
#
#
#r=0.2
#
#rect= Rectangle(Point(xinf,yinf),Point(xsup,ysup))
#circle = Circle(Point(c_x,c_y),0.45)
#circle1 = Circle(Point(c1_x,c1_y),r)
#domain=rect-(circle+circle1)
#res = 50  # Resolution of mesh
#mesh = generate_mesh(domain, res)
#
##### Premier raffinement dans l axe horizontale du cyٍlindre 
#meshB= raffinemment_maillage_q_1(r,mesh)
#mesh=meshB
#
#class Obstacle(SubDomain):
#    def inside(self, x, on_boundary):
#        return (on_boundary and between((x[0]-c_x), (-0.45-tol, 0.45+tol)) and between((x[1]-c_x), (-0.45-tol, 0.45+tol)))
#
#
#class Obstacle1(SubDomain):
#    def inside(self, x, on_boundary):
#        return (on_boundary and between((x[0]-c1_x), (-r-tol, r+tol)) and between((x[1]-c1_y), (-r-tol, r+tol)))
#
#
#
#V=VectorFunctionSpace(mesh, 'P', 2, constrained_domain=PeriodicBoundary())
#nb_noeuds_ess=V.dim()
#
#print ('uuuuuuu',nb_noeuds_ess)
#
#obstacle = Obstacle()
#obstacle1 = Obstacle1()
#
#boundaries = FacetFunction("size_t", mesh)
#
#boundaries.set_all(5)
#obstacle.mark(boundaries, 0)
#obstacle1.mark(boundaries, 1)
#
#    
#ds = Measure("ds")[boundaries]
#normale = FacetNormal(mesh) ## normale associe a chaque cellule
#
#num0=0
#num1=1
#
##################################################
######Sauvegarder les phi apres interpolation#####
##################################################  
#w=np.zeros((nb_noeuds_ess,nb_modes))
#wi_fixe=Function(V_fixe)
#
#repertoire_final_1="Homog_base_maill_init_multi/"
#
#File("Homog_base_maill_init_multi/mesh_quelconque.xml.gz") << mesh
#
#ufile,USAVE=creation_fichier_pourecriture_champ_hdf5(repertoire_final_1,mesh)
#file_rayon_ecriture = open("%s/rayon_ecriture.txt" %(repertoire_final_1), "w")
#kfic=0
#for n in range(nb_modes):
#    xc=n+1
#    wi_fixe.vector().set_local(G[:,n])
#    wi_fixe.set_allow_extrapolation(True)
#    wi_ess=interpolate(wi_fixe,V)
##    plot(wi_ess)
##    plt.show()
#    ecriture_champ_hdf5(ufile,USAVE,wi_ess,kfic,file_rayon_ecriture,xc)
#    kfic=kfic+1
#    w[:,n]=wi_ess.vector().array()
#    
#
#
######Model reduit pour un rayon précis#####
#mat_mass=np.zeros((nb_modes,nb_modes))
#vec_sf=np.zeros((nb_modes))
#phiu_i = Function(V)                                                                    
#phiu_j = Function(V)                                                                    
#    
#
#
#for iproj in range(nb_modes):
#    phiu_i.vector().set_local(w[:,iproj])
#	
#        
#    for jrec in range(nb_modes):
#          phiu_j.vector().set_local(w[:,jrec])
#              
#              
#          ### matrice de masse
#          mat_mass[jrec,iproj]=assemble(tr(dot((grad(phiu_j)).T, grad(phiu_i)))*dx)
#          
#
#for iproj in range(nb_modes):
#    phiu_i.vector().set_local(w[:,iproj])
#    ### terme sf
#    vec_sf[iproj]=assemble((dot(normale,phiu_i)*ds(num0))+(dot(normale,phiu_i)*ds(num1)))
#
#sol=np.linalg.solve(mat_mass.T,-vec_sf)
#
#print sol
#
#repertoire_final_2="Solution_homog_red_maill_init_multi/"
#
#ufile,USAVE=creation_fichier_pourecriture_champ_hdf5(repertoire_final_2,mesh)
#file_rayon_ecriture = open("%s/rayon_ecriture.txt" %(repertoire_final_2), "w")
#kfic=0
#
#
#urec=np.dot(w,sol)
#urec_red=Function(V)
#urec_red.vector().set_local(urec)
#
#ecriture_champ_hdf5(ufile,USAVE,urec_red,kfic,file_rayon_ecriture,1)
#
#H=assemble(grad(urec_red)[0,0]*dx)
#M=assemble(grad(urec_red)[1,1]*dx)
#X=assemble(grad(urec_red)[0,1]*dx)
#C=assemble(grad(urec_red)[1,0]*dx)
#print H,M,X,C
#ep=(1-pi*(0.25*(r**2)+0.45**2))
#D=((H/ep)+1)*ep
#print ("hi",r,ep,D)
#
#
#
#plot(urec_red)
#plt.show()





#################################################################
################### 3 quart et cercle central####################
#################################################################
#########Interpolation de la base dans chaque maillage###########
#################################################################
#
#
#r=0.2
#
#rect= Rectangle(Point(xinf,yinf),Point(xsup,ysup))
#circle = Circle(Point(c_x,c_y),0.45)
#circle1 = Circle(Point(c1_x,c1_y),r)
#circle2 = Circle(Point(c2_x,c2_y),r)
#circle3 = Circle(Point(c3_x,c3_y),r)
#domain=rect-(circle+circle1+circle2+circle3)
#res = 50  # Resolution of mesh
#mesh = generate_mesh(domain, res)
#
##### Premier raffinement dans l axe horizontale du cyٍlindre 
#meshB= raffinemment_maillage_q_2(r,mesh)
#mesh=meshB
#
#class Obstacle(SubDomain):
#    def inside(self, x, on_boundary):
#        return (on_boundary and between((x[0]-c_x), (-0.45-tol, 0.45+tol)) and between((x[1]-c_x), (-0.45-tol, 0.45+tol)))
#
#
#class Obstacle1(SubDomain):
#    def inside(self, x, on_boundary):
#        return (on_boundary and between((x[0]-c1_x), (-r-tol, r+tol)) and between((x[1]-c1_y), (-r-tol, r+tol)))
#
#
#class Obstacle2(SubDomain):
#    def inside(self, x, on_boundary):
#        return (on_boundary and between((x[0]-c2_x), (-r-tol, r+tol)) and between((x[1]-c2_y), (-r-tol, r+tol)))
#
#
#class Obstacle3(SubDomain):
#    def inside(self, x, on_boundary):
#        return (on_boundary and between((x[0]-c3_x), (-r-tol, r+tol)) and between((x[1]-c3_y), (-r-tol, r+tol)))
#
#
#V=VectorFunctionSpace(mesh, 'P', 2, constrained_domain=PeriodicBoundary())
#nb_noeuds_ess=V.dim()
#
#print ('uuuuuuu',nb_noeuds_ess)
#
#obstacle = Obstacle()
#obstacle1 = Obstacle1()
#obstacle2 = Obstacle2()
#obstacle3 = Obstacle3()
#
#
#boundaries = FacetFunction("size_t", mesh)
#
#boundaries.set_all(5)
#obstacle.mark(boundaries, 0)
#obstacle1.mark(boundaries, 1)
#obstacle2.mark(boundaries, 2)
#obstacle3.mark(boundaries, 3)
#
#    
#ds = Measure("ds")[boundaries]
#normale = FacetNormal(mesh) ## normale associe a chaque cellule
#
#num0=0
#num1=1
#num2=2
#num3=3
#
##################################################
######Sauvegarder les phi apres interpolation#####
##################################################  
#w=np.zeros((nb_noeuds_ess,nb_modes))
#wi_fixe=Function(V_fixe)
#
#repertoire_final_1="Homog_base_maill_init_multi/"
#
#File("Homog_base_maill_init_multi/mesh_quelconque.xml.gz") << mesh
#
#ufile,USAVE=creation_fichier_pourecriture_champ_hdf5(repertoire_final_1,mesh)
#file_rayon_ecriture = open("%s/rayon_ecriture.txt" %(repertoire_final_1), "w")
#kfic=0
#for n in range(nb_modes):
#    xc=n+1
#    wi_fixe.vector().set_local(G[:,n])
#    wi_fixe.set_allow_extrapolation(True)
#    wi_ess=interpolate(wi_fixe,V)
##    plot(wi_ess)
##    plt.show()
#    ecriture_champ_hdf5(ufile,USAVE,wi_ess,kfic,file_rayon_ecriture,xc)
#    kfic=kfic+1
#    w[:,n]=wi_ess.vector().array()
#    
#
#
######Model reduit pour un rayon précis#####
#mat_mass=np.zeros((nb_modes,nb_modes))
#vec_sf=np.zeros((nb_modes))
#phiu_i = Function(V)                                                                    
#phiu_j = Function(V)                                                                    
#    
#
#
#for iproj in range(nb_modes):
#    phiu_i.vector().set_local(w[:,iproj])
#	
#        
#    for jrec in range(nb_modes):
#          phiu_j.vector().set_local(w[:,jrec])
#              
#              
#          ### matrice de masse
#          mat_mass[jrec,iproj]=assemble(tr(dot((grad(phiu_j)).T, grad(phiu_i)))*dx)
#          
#
#for iproj in range(nb_modes):
#    phiu_i.vector().set_local(w[:,iproj])
#    ### terme sf
#    vec_sf[iproj]=assemble((dot(normale,phiu_i)*ds(num0))+(dot(normale,phiu_i)*ds(num1))+(dot(normale,phiu_i)*ds(num2))+(dot(normale,phiu_i)*ds(num3)))
#
#sol=np.linalg.solve(mat_mass.T,-vec_sf)
#
#print sol
#
#repertoire_final_2="Solution_homog_red_maill_init_multi/"
#
#ufile,USAVE=creation_fichier_pourecriture_champ_hdf5(repertoire_final_2,mesh)
#file_rayon_ecriture = open("%s/rayon_ecriture.txt" %(repertoire_final_2), "w")
#kfic=0
#
#
#urec=np.dot(w,sol)
#urec_red=Function(V)
#urec_red.vector().set_local(urec)
#
#ecriture_champ_hdf5(ufile,USAVE,urec_red,kfic,file_rayon_ecriture,1)
#
#H=assemble(grad(urec_red)[0,0]*dx)
#M=assemble(grad(urec_red)[1,1]*dx)
#X=assemble(grad(urec_red)[0,1]*dx)
#C=assemble(grad(urec_red)[1,0]*dx)
#print H,M,X,C
#ep=(1-pi*(0.75*(r**2)+0.45**2))
#D=((H/ep)+1)*ep
#print ("hi",r,ep,D)
#
#
#
#plot(urec_red)
#plt.show()








#################################################################
######### 2 quart sur même axe et cercle central#################
#################################################################
#########Interpolation de la base dans chaque maillage###########
#################################################################
#
#
#r=0.2
#
#rect= Rectangle(Point(xinf,yinf),Point(xsup,ysup))
#circle = Circle(Point(c_x,c_y),0.45)
#circle1 = Circle(Point(c1_x,c1_y),r)
#circle2 = Circle(Point(c2_x,c2_y),r)
#domain=rect-(circle+circle1+circle2)
#res = 50  # Resolution of mesh
#mesh = generate_mesh(domain, res)
#
##### Premier raffinement dans l axe horizontale du cyٍlindre 
#meshB= raffinemment_maillage_q_3(r,mesh)
#mesh=meshB
#
#class Obstacle(SubDomain):
#    def inside(self, x, on_boundary):
#        return (on_boundary and between((x[0]-c_x), (-0.45-tol, 0.45+tol)) and between((x[1]-c_x), (-0.45-tol, 0.45+tol)))
#
#
#class Obstacle1(SubDomain):
#    def inside(self, x, on_boundary):
#        return (on_boundary and between((x[0]-c1_x), (-r-tol, r+tol)) and between((x[1]-c1_y), (-r-tol, r+tol)))
#
#
#class Obstacle2(SubDomain):
#    def inside(self, x, on_boundary):
#        return (on_boundary and between((x[0]-c2_x), (-r-tol, r+tol)) and between((x[1]-c2_y), (-r-tol, r+tol)))
#
#
#
#
#V=VectorFunctionSpace(mesh, 'P', 2, constrained_domain=PeriodicBoundary())
#nb_noeuds_ess=V.dim()
#
#print ('uuuuuuu',nb_noeuds_ess)
#
#obstacle = Obstacle()
#obstacle1 = Obstacle1()
#obstacle2 = Obstacle2()
#
#
#boundaries = FacetFunction("size_t", mesh)
#
#boundaries.set_all(5)
#obstacle.mark(boundaries, 0)
#obstacle1.mark(boundaries, 1)
#obstacle2.mark(boundaries, 2)
#
#    
#ds = Measure("ds")[boundaries]
#normale = FacetNormal(mesh) ## normale associe a chaque cellule
#
#num0=0
#num1=1
#num2=2
#
##################################################
######Sauvegarder les phi apres interpolation#####
##################################################  
#w=np.zeros((nb_noeuds_ess,nb_modes))
#wi_fixe=Function(V_fixe)
#
#repertoire_final_1="Homog_base_maill_init_multi/"
#
#File("Homog_base_maill_init_multi/mesh_quelconque.xml.gz") << mesh
#
#ufile,USAVE=creation_fichier_pourecriture_champ_hdf5(repertoire_final_1,mesh)
#file_rayon_ecriture = open("%s/rayon_ecriture.txt" %(repertoire_final_1), "w")
#kfic=0
#for n in range(nb_modes):
#    xc=n+1
#    wi_fixe.vector().set_local(G[:,n])
#    wi_fixe.set_allow_extrapolation(True)
#    wi_ess=interpolate(wi_fixe,V)
##    plot(wi_ess)
##    plt.show()
#    ecriture_champ_hdf5(ufile,USAVE,wi_ess,kfic,file_rayon_ecriture,xc)
#    kfic=kfic+1
#    w[:,n]=wi_ess.vector().array()
#    
#
#
######Model reduit pour un rayon précis#####
#mat_mass=np.zeros((nb_modes,nb_modes))
#vec_sf=np.zeros((nb_modes))
#phiu_i = Function(V)                                                                    
#phiu_j = Function(V)                                                                    
#    
#
#
#for iproj in range(nb_modes):
#    phiu_i.vector().set_local(w[:,iproj])
#	
#        
#    for jrec in range(nb_modes):
#          phiu_j.vector().set_local(w[:,jrec])
#              
#              
#          ### matrice de masse
#          mat_mass[jrec,iproj]=assemble(tr(dot((grad(phiu_j)).T, grad(phiu_i)))*dx)
#          
#
#for iproj in range(nb_modes):
#    phiu_i.vector().set_local(w[:,iproj])
#    ### terme sf
#    vec_sf[iproj]=assemble((dot(normale,phiu_i)*ds(num0))+(dot(normale,phiu_i)*ds(num1))+(dot(normale,phiu_i)*ds(num2)))
#
#sol=np.linalg.solve(mat_mass.T,-vec_sf)
#
#print sol
#
#repertoire_final_2="Solution_homog_red_maill_init_multi/"
#
#ufile,USAVE=creation_fichier_pourecriture_champ_hdf5(repertoire_final_2,mesh)
#file_rayon_ecriture = open("%s/rayon_ecriture.txt" %(repertoire_final_2), "w")
#kfic=0
#
#
#urec=np.dot(w,sol)
#urec_red=Function(V)
#urec_red.vector().set_local(urec)
#
#ecriture_champ_hdf5(ufile,USAVE,urec_red,kfic,file_rayon_ecriture,1)
#
#H=assemble(grad(urec_red)[0,0]*dx)
#M=assemble(grad(urec_red)[1,1]*dx)
#X=assemble(grad(urec_red)[0,1]*dx)
#C=assemble(grad(urec_red)[1,0]*dx)
#print H,M,X,C
#ep=(1-pi*(0.5*(r**2)+0.45**2))
#D=((H/ep)+1)*ep
#print ("hi",r,ep,D)
#
#
#
#plot(urec_red)
#plt.show()








#################################################################
######### 2 quart sur diagonal et cercle central#################
#################################################################
#########Interpolation de la base dans chaque maillage###########
#################################################################
#
#
#r=0.2
#
#rect= Rectangle(Point(xinf,yinf),Point(xsup,ysup))
#circle = Circle(Point(c_x,c_y),0.45)
#circle1 = Circle(Point(c1_x,c1_y),r)
#circle3 = Circle(Point(c3_x,c3_y),r)
#domain=rect-(circle+circle1+circle3)
#res = 50  # Resolution of mesh
#mesh = generate_mesh(domain, res)
#
##### Premier raffinement dans l axe horizontale du cyٍlindre 
#meshB= raffinemment_maillage_q_4(r,mesh)
#mesh=meshB
#
#class Obstacle(SubDomain):
#    def inside(self, x, on_boundary):
#        return (on_boundary and between((x[0]-c_x), (-0.45-tol, 0.45+tol)) and between((x[1]-c_x), (-0.45-tol, 0.45+tol)))
#
#
#class Obstacle1(SubDomain):
#    def inside(self, x, on_boundary):
#        return (on_boundary and between((x[0]-c1_x), (-r-tol, r+tol)) and between((x[1]-c1_y), (-r-tol, r+tol)))
#
#
#class Obstacle3(SubDomain):
#    def inside(self, x, on_boundary):
#        return (on_boundary and between((x[0]-c3_x), (-r-tol, r+tol)) and between((x[1]-c3_y), (-r-tol, r+tol)))
#
#
#
#
#V=VectorFunctionSpace(mesh, 'P', 2, constrained_domain=PeriodicBoundary())
#nb_noeuds_ess=V.dim()
#
#print ('uuuuuuu',nb_noeuds_ess)
#
#obstacle = Obstacle()
#obstacle1 = Obstacle1()
#obstacle3 = Obstacle3()
#
#
#boundaries = FacetFunction("size_t", mesh)
#
#boundaries.set_all(5)
#obstacle.mark(boundaries, 0)
#obstacle1.mark(boundaries, 1)
#obstacle3.mark(boundaries, 3)
#
#    
#ds = Measure("ds")[boundaries]
#normale = FacetNormal(mesh) ## normale associe a chaque cellule
#
#num0=0
#num1=1
#num3=3
#
##################################################
######Sauvegarder les phi apres interpolation#####
##################################################  
#w=np.zeros((nb_noeuds_ess,nb_modes))
#wi_fixe=Function(V_fixe)
#
#repertoire_final_1="Homog_base_maill_init_multi/"
#
#File("Homog_base_maill_init_multi/mesh_quelconque.xml.gz") << mesh
#
#ufile,USAVE=creation_fichier_pourecriture_champ_hdf5(repertoire_final_1,mesh)
#file_rayon_ecriture = open("%s/rayon_ecriture.txt" %(repertoire_final_1), "w")
#kfic=0
#for n in range(nb_modes):
#    xc=n+1
#    wi_fixe.vector().set_local(G[:,n])
#    wi_fixe.set_allow_extrapolation(True)
#    wi_ess=interpolate(wi_fixe,V)
##    plot(wi_ess)
##    plt.show()
#    ecriture_champ_hdf5(ufile,USAVE,wi_ess,kfic,file_rayon_ecriture,xc)
#    kfic=kfic+1
#    w[:,n]=wi_ess.vector().array()
#    
#
#
######Model reduit pour un rayon précis#####
#mat_mass=np.zeros((nb_modes,nb_modes))
#vec_sf=np.zeros((nb_modes))
#phiu_i = Function(V)                                                                    
#phiu_j = Function(V)                                                                    
#    
#
#
#for iproj in range(nb_modes):
#    phiu_i.vector().set_local(w[:,iproj])
#	
#        
#    for jrec in range(nb_modes):
#          phiu_j.vector().set_local(w[:,jrec])
#              
#              
#          ### matrice de masse
#          mat_mass[jrec,iproj]=assemble(tr(dot((grad(phiu_j)).T, grad(phiu_i)))*dx)
#          
#
#for iproj in range(nb_modes):
#    phiu_i.vector().set_local(w[:,iproj])
#    ### terme sf
#    vec_sf[iproj]=assemble((dot(normale,phiu_i)*ds(num0))+(dot(normale,phiu_i)*ds(num1))+(dot(normale,phiu_i)*ds(num3)))
#
#sol=np.linalg.solve(mat_mass.T,-vec_sf)
#
#print sol
#
#repertoire_final_2="Solution_homog_red_maill_init_multi/"
#
#ufile,USAVE=creation_fichier_pourecriture_champ_hdf5(repertoire_final_2,mesh)
#file_rayon_ecriture = open("%s/rayon_ecriture.txt" %(repertoire_final_2), "w")
#kfic=0
#
#
#urec=np.dot(w,sol)
#urec_red=Function(V)
#urec_red.vector().set_local(urec)
#
#ecriture_champ_hdf5(ufile,USAVE,urec_red,kfic,file_rayon_ecriture,1)
#
#H=assemble(grad(urec_red)[0,0]*dx)
#M=assemble(grad(urec_red)[1,1]*dx)
#X=assemble(grad(urec_red)[0,1]*dx)
#C=assemble(grad(urec_red)[1,0]*dx)
#print H,M,X,C
#ep=(1-pi*(0.5*(r**2)+0.45**2))
#D=((H/ep)+1)*ep
#print ("hi",r,ep,D)
#
#
#
#plot(urec_red)
#plt.show()










#################################################################
###########4 differents quart et cercle central##################
#################################################################
#########Interpolation de la base dans chaque maillage###########
#################################################################
#r1=0.05
#r2=0.1
#r3=0.15
#r4=0.2
##r1=r2=r3=r4
#
#
#rect= Rectangle(Point(xinf,yinf),Point(xsup,ysup))
#circle = Circle(Point(c_x,c_y),0.45)
#circle1 = Circle(Point(c1_x,c1_y),r1)
#circle2 = Circle(Point(c2_x,c2_y),r2)
#circle3 = Circle(Point(c3_x,c3_y),r3)
#circle4 = Circle(Point(c4_x,c4_y),r4)
#domain=rect-(circle+circle1+circle2+circle3+circle4)
#res = 50  # Resolution of mesh
#mesh = generate_mesh(domain, res)
#
##### Premier raffinement dans l axe horizontale du cyٍlindre 
#meshB= raffinemment_maillage_q_5(r1,r2,r3,r4,mesh)
#mesh=meshB
#
#class Obstacle(SubDomain):
#    def inside(self, x, on_boundary):
#        return (on_boundary and between((x[0]-c_x), (-0.45-tol, 0.45+tol)) and between((x[1]-c_x), (-0.45-tol, 0.45+tol)))
#
#
#class Obstacle1(SubDomain):
#    def inside(self, x, on_boundary):
#        return (on_boundary and between((x[0]-c1_x), (-r1-tol, r1+tol)) and between((x[1]-c1_y), (-r1-tol, r1+tol)))
#
#
#class Obstacle2(SubDomain):
#    def inside(self, x, on_boundary):
#        return (on_boundary and between((x[0]-c2_x), (-r2-tol, r2+tol)) and between((x[1]-c2_y), (-r2-tol, r2+tol)))
#
#
#class Obstacle3(SubDomain):
#    def inside(self, x, on_boundary):
#        return (on_boundary and between((x[0]-c3_x), (-r3-tol, r3+tol)) and between((x[1]-c3_y), (-r3-tol, r3+tol)))
#
#
#class Obstacle4(SubDomain):
#    def inside(self, x, on_boundary):
#        return (on_boundary and between((x[0]-c4_x), (-r4-tol, r4+tol)) and between((x[1]-c4_y), (-r4-tol, r4+tol)))
#
#
#V=VectorFunctionSpace(mesh, 'P', 2, constrained_domain=PeriodicBoundary())
#nb_noeuds_ess=V.dim()
#
#print ('uuuuuuu',nb_noeuds_ess)
#
#obstacle = Obstacle()
#obstacle1 = Obstacle1()
#obstacle2 = Obstacle2()
#obstacle3 = Obstacle3()
#obstacle4 = Obstacle4()
#
#boundaries = FacetFunction("size_t", mesh)
#
#boundaries.set_all(5)
#obstacle.mark(boundaries, 0)
#obstacle1.mark(boundaries, 1)
#obstacle2.mark(boundaries, 2)
#obstacle3.mark(boundaries, 3)
#obstacle4.mark(boundaries, 4)
#
#    
#ds = Measure("ds")[boundaries]
#normale = FacetNormal(mesh) ## normale associe a chaque cellule
#
#num0=0
#num1=1
#num2=2
#num3=3
#num4=4
#
#
##################################################
######Sauvegarder les phi apres interpolation#####
##################################################  
#w=np.zeros((nb_noeuds_ess,nb_modes))
#wi_fixe=Function(V_fixe)
#
#repertoire_final_1="Homog_base_maill_init_multi/"
#
#File("Homog_base_maill_init_multi/mesh_complexe.xml.gz") << mesh
#
#ufile,USAVE=creation_fichier_pourecriture_champ_hdf5(repertoire_final_1,mesh)
#file_rayon_ecriture = open("%s/rayon_ecriture.txt" %(repertoire_final_1), "w")
#kfic=0
#for n in range(nb_modes):
#    xc=n+1
#    wi_fixe.vector().set_local(G[:,n])
#    wi_fixe.set_allow_extrapolation(True)
#    wi_ess=interpolate(wi_fixe,V)
##    plot(wi_ess)
##    plt.show()
#    ecriture_champ_hdf5(ufile,USAVE,wi_ess,kfic,file_rayon_ecriture,xc)
#    kfic=kfic+1
#    w[:,n]=wi_ess.vector().array()
#    
#
#
######Model reduit pour un rayon précis#####
#mat_mass=np.zeros((nb_modes,nb_modes))
#vec_sf=np.zeros((nb_modes))
#phiu_i = Function(V)                                                                    
#phiu_j = Function(V)                                                                    
#    
#
#
#for iproj in range(nb_modes):
#    phiu_i.vector().set_local(w[:,iproj])
#	
#        
#    for jrec in range(nb_modes):
#          phiu_j.vector().set_local(w[:,jrec])
#              
#              
#          ### matrice de masse
#          mat_mass[jrec,iproj]=assemble(tr(dot((grad(phiu_j)).T, grad(phiu_i)))*dx)
#          
#
#for iproj in range(nb_modes):
#    phiu_i.vector().set_local(w[:,iproj])
#    ### terme sf
#    vec_sf[iproj]=assemble((dot(normale,phiu_i)*ds(num0))+(dot(normale,phiu_i)*ds(num1))+(dot(normale,phiu_i)*ds(num2))+(dot(normale,phiu_i)*ds(num3))+(dot(normale,phiu_i)*ds(num4)))
#
#sol=np.linalg.solve(mat_mass.T,-vec_sf)
#
#print sol
#
#repertoire_final_2="Solution_homog_red_maill_init_multi/"
#
#ufile,USAVE=creation_fichier_pourecriture_champ_hdf5(repertoire_final_2,mesh)
#file_rayon_ecriture = open("%s/rayon_ecriture.txt" %(repertoire_final_2), "w")
#kfic=0
#
#
#urec=np.dot(w,sol)
#urec_red=Function(V)
#urec_red.vector().set_local(urec)
#
#ecriture_champ_hdf5(ufile,USAVE,urec_red,kfic,file_rayon_ecriture,1)
#
#H=assemble(grad(urec_red)[0,0]*dx)
#M=assemble(grad(urec_red)[1,1]*dx)
#X=assemble(grad(urec_red)[0,1]*dx)
#C=assemble(grad(urec_red)[1,0]*dx)
#print H,M,X,C
#ep=(1-pi*(0.25*(r1**2)+0.25*(r2**2)+0.25*(r3**2)+0.25*(r4**2)+0.45**2))
#D=((H/ep)+1)*ep
#print ("hi",r,ep,D)
#
#
#
#plot(urec_red)
#plt.show()











#################################################################
######################## 4 différentes quarts####################
#################################################################
#########Interpolation de la base dans chaque maillage###########
#################################################################
#r1=0.1
#r2=0.2
#r3=0.3
#r4=0.4
#
#rect= Rectangle(Point(xinf,yinf),Point(xsup,ysup))
#circle1 = Circle(Point(c1_x,c1_y),r1)
#circle2 = Circle(Point(c2_x,c2_y),r2)
#circle3 = Circle(Point(c3_x,c3_y),r3)
#circle4 = Circle(Point(c4_x,c4_y),r4)
#domain=rect-(circle1+circle2+circle3+circle4)
#res = 50  # Resolution of mesh
#mesh = generate_mesh(domain, res)
#
##### Premier raffinement dans l axe horizontale du cyٍlindre 
#meshB= raffinemment_maillage_q_6(r1,r2,r3,r4,mesh)
#mesh=meshB
#
#
#class Obstacle1(SubDomain):
#    def inside(self, x, on_boundary):
#        return (on_boundary and between((x[0]-c1_x), (-r1-tol, r1+tol)) and between((x[1]-c1_y), (-r1-tol, r1+tol)))
#
#
#class Obstacle2(SubDomain):
#    def inside(self, x, on_boundary):
#        return (on_boundary and between((x[0]-c2_x), (-r2-tol, r2+tol)) and between((x[1]-c2_y), (-r2-tol, r2+tol)))
#
#
#class Obstacle3(SubDomain):
#    def inside(self, x, on_boundary):
#        return (on_boundary and between((x[0]-c3_x), (-r3-tol, r3+tol)) and between((x[1]-c3_y), (-r3-tol, r3+tol)))
#
#
#class Obstacle4(SubDomain):
#    def inside(self, x, on_boundary):
#        return (on_boundary and between((x[0]-c4_x), (-r4-tol, r4+tol)) and between((x[1]-c4_y), (-r4-tol, r4+tol)))
#
#
#
#V=VectorFunctionSpace(mesh, 'P', 2, constrained_domain=PeriodicBoundary())
#nb_noeuds_ess=V.dim()
#
#print ('uuuuuuu',nb_noeuds_ess)
#
#obstacle1 = Obstacle1()
#obstacle2 = Obstacle2()
#obstacle3 = Obstacle3()
#obstacle4 = Obstacle4()
#
#boundaries = FacetFunction("size_t", mesh)
#
#boundaries.set_all(0)
#
#obstacle1.mark(boundaries, 1)
#obstacle2.mark(boundaries, 2)
#obstacle3.mark(boundaries, 3)
#obstacle4.mark(boundaries, 4)
#
#    
#ds = Measure("ds")[boundaries]
#normale = FacetNormal(mesh) ## normale associe a chaque cellule
#
#
#num1=1
#num2=2
#num3=3
#num4=4
#
#
##################################################
######Sauvegarder les phi apres interpolation#####
##################################################  
#w=np.zeros((nb_noeuds_ess,nb_modes))
#wi_fixe=Function(V_fixe)
#
#repertoire_final_1="Homog_base_maill_init_multi/"
#
#File("Homog_base_maill_init_multi/mesh_complexe.xml.gz") << mesh
#
#ufile,USAVE=creation_fichier_pourecriture_champ_hdf5(repertoire_final_1,mesh)
#file_rayon_ecriture = open("%s/rayon_ecriture.txt" %(repertoire_final_1), "w")
#kfic=0
#for n in range(nb_modes):
#    xc=n+1
#    wi_fixe.vector().set_local(G[:,n])
#    wi_fixe.set_allow_extrapolation(True)
#    wi_ess=interpolate(wi_fixe,V)
##    plot(wi_ess)
##    plt.show()
#    ecriture_champ_hdf5(ufile,USAVE,wi_ess,kfic,file_rayon_ecriture,xc)
#    kfic=kfic+1
#    w[:,n]=wi_ess.vector().array()
#    
#
#
######Model reduit pour un rayon précis#####
#mat_mass=np.zeros((nb_modes,nb_modes))
#vec_sf=np.zeros((nb_modes))
#phiu_i = Function(V)                                                                    
#phiu_j = Function(V)                                                                    
#    
#
#
#for iproj in range(nb_modes):
#    phiu_i.vector().set_local(w[:,iproj])
#	
#        
#    for jrec in range(nb_modes):
#          phiu_j.vector().set_local(w[:,jrec])
#              
#              
#          ### matrice de masse
#          mat_mass[jrec,iproj]=assemble(tr(dot((grad(phiu_j)).T, grad(phiu_i)))*dx)
#          
#
#for iproj in range(nb_modes):
#    phiu_i.vector().set_local(w[:,iproj])
#    ### terme sf
#    vec_sf[iproj]=assemble((dot(normale,phiu_i)*ds(num1))+(dot(normale,phiu_i)*ds(num2))+(dot(normale,phiu_i)*ds(num3))+(dot(normale,phiu_i)*ds(num4)))
#
#sol=np.linalg.solve(mat_mass.T,-vec_sf)
#
#print sol
#
#repertoire_final_2="Solution_homog_red_maill_init_multi/"
#
#ufile,USAVE=creation_fichier_pourecriture_champ_hdf5(repertoire_final_2,mesh)
#file_rayon_ecriture = open("%s/rayon_ecriture.txt" %(repertoire_final_2), "w")
#kfic=0
#
#
#urec=np.dot(w,sol)
#urec_red=Function(V)
#urec_red.vector().set_local(urec)
#
#ecriture_champ_hdf5(ufile,USAVE,urec_red,kfic,file_rayon_ecriture,1)
#
#H=assemble(grad(urec_red)[0,0]*dx)
#M=assemble(grad(urec_red)[1,1]*dx)
#X=assemble(grad(urec_red)[0,1]*dx)
#C=assemble(grad(urec_red)[1,0]*dx)
#print H,M,X,C
#ep=(1-pi*(0.25*r1**2+0.25*r2**2+0.25*r3**2+0.25*r4**2))
#D=((H/ep)+1)*ep
#print ("hi",r,ep,D)
#
##print urec_red(0,0.5),urec_red(1,0.5),urec_red(0.5,0),urec_red(0.5,1),urec_red(0,0.4),urec_red(1,0.4)
#
#plot(urec_red)
#plt.show()











################################################################
################ 4 rectangle avec cercle########################
################################################################
########Interpolation de la base dans chaque maillage###########
################################################################



rect= Rectangle(Point(xinf,yinf),Point(xsup,ysup))
circle = Circle(Point(c_x,c_y),0.45)
circle1 = Rectangle(Point(xinf,yinf),Point(0.12,0.12))
circle2 = Rectangle(Point(0.88,0.12),Point(1.0,0))
circle3 = Rectangle(Point(0.88,0.88),Point(xsup,ysup))
circle4 = Rectangle(Point(0.12,0.88),Point(0,1.0))
domain=rect-(circle+circle1+circle2+circle3+circle4)
res = 80  # Resolution of mesh
mesh = generate_mesh(domain, res)

plot(mesh)
plt.show()

class Obstacle(SubDomain):
    def inside(self, x, on_boundary):
        return (on_boundary and between((x[0]-c_x), (-0.45-tol, 0.45+tol)) and between((x[1]-c_x), (-0.45-tol, 0.45+tol)))


class Obstacle1(SubDomain):
    def inside(self, x, on_boundary):
        return (on_boundary and between(x[0], (-0.1,0.13)) and between(x[1], (-0.1,0.13)))


class Obstacle2(SubDomain):
    def inside(self, x, on_boundary):
        return (on_boundary and between(x[0], (0.87,1.1)) and between(x[1], (-0.1,0.13)))


class Obstacle3(SubDomain):
    def inside(self, x, on_boundary):
        return (on_boundary and between(x[0], (0.87,1.1)) and between(x[1], (0.87,1.1)))


class Obstacle4(SubDomain):
    def inside(self, x, on_boundary):
        return (on_boundary and between(x[0], (-0.1,0.13)) and between(x[1], (0.87,1.1)))


V=VectorFunctionSpace(mesh, 'P', 2, constrained_domain=PeriodicBoundary())
nb_noeuds_ess=V.dim()

print ('uuuuuuu',nb_noeuds_ess)

obstacle = Obstacle()
obstacle1 = Obstacle1()
obstacle2 = Obstacle2()
obstacle3 = Obstacle3()
obstacle4 = Obstacle4()

boundaries = FacetFunction("size_t", mesh)

boundaries.set_all(5)
obstacle.mark(boundaries, 0)
obstacle1.mark(boundaries, 1)
obstacle2.mark(boundaries, 2)
obstacle3.mark(boundaries, 3)
obstacle4.mark(boundaries, 4)

    
ds = Measure("ds")[boundaries]
normale = FacetNormal(mesh) ## normale associe a chaque cellule

num0=0
num1=1
num2=2
num3=3
num4=4


#################################################
#####Sauvegarder les phi apres interpolation#####
#################################################  
w=np.zeros((nb_noeuds_ess,nb_modes))
wi_fixe=Function(V_fixe)

repertoire_final_1="Homog_base_maill_init_multi/"

File("Homog_base_maill_init_multi/mesh_complexe.xml.gz") << mesh

ufile,USAVE=creation_fichier_pourecriture_champ_hdf5(repertoire_final_1,mesh)
file_rayon_ecriture = open("%s/rayon_ecriture.txt" %(repertoire_final_1), "w")
kfic=0
for n in range(nb_modes):
    xc=n+1
    wi_fixe.vector().set_local(G[:,n])
    wi_fixe.set_allow_extrapolation(True)
    wi_ess=interpolate(wi_fixe,V)
#    plot(wi_ess)
#    plt.show()
    ecriture_champ_hdf5(ufile,USAVE,wi_ess,kfic,file_rayon_ecriture,xc)
    kfic=kfic+1
    w[:,n]=wi_ess.vector().array()
    


#####Model reduit pour un rayon précis#####
mat_mass=np.zeros((nb_modes,nb_modes))
vec_sf=np.zeros((nb_modes))
phiu_i = Function(V)                                                                    
phiu_j = Function(V)                                                                    
    


for iproj in range(nb_modes):
    phiu_i.vector().set_local(w[:,iproj])
	
        
    for jrec in range(nb_modes):
          phiu_j.vector().set_local(w[:,jrec])
              
              
          ### matrice de masse
          mat_mass[jrec,iproj]=assemble(tr(dot((grad(phiu_j)).T, grad(phiu_i)))*dx)
          

for iproj in range(nb_modes):
    phiu_i.vector().set_local(w[:,iproj])
    ### terme sf
    vec_sf[iproj]=assemble((dot(normale,phiu_i)*ds(num0))+(dot(normale,phiu_i)*ds(num1))+(dot(normale,phiu_i)*ds(num2))+(dot(normale,phiu_i)*ds(num3))+(dot(normale,phiu_i)*ds(num4)))

sol=np.linalg.solve(mat_mass.T,-vec_sf)

print sol

repertoire_final_2="Solution_homog_red_maill_init_multi/"

ufile,USAVE=creation_fichier_pourecriture_champ_hdf5(repertoire_final_2,mesh)
file_rayon_ecriture = open("%s/rayon_ecriture.txt" %(repertoire_final_2), "w")
kfic=0


urec=np.dot(w,sol)
urec_red=Function(V)
urec_red.vector().set_local(urec)

ecriture_champ_hdf5(ufile,USAVE,urec_red,kfic,file_rayon_ecriture,1)

H=assemble(grad(urec_red)[0,0]*dx)
M=assemble(grad(urec_red)[1,1]*dx)
X=assemble(grad(urec_red)[0,1]*dx)
C=assemble(grad(urec_red)[1,0]*dx)
print H,M,X,C
ep=(1-pi*(0.45**2)-4*0.12*0.12)
D=((H/ep)+1)*ep
print ("hi",r,ep,D)



plot(urec_red)
plt.show()










#################################################################
################# 4 rectangle avec rectangle central#############
#################################################################
#########Interpolation de la base dans chaque maillage###########
#################################################################
#
#
#
#rect= Rectangle(Point(xinf,yinf),Point(xsup,ysup))
#circle = Rectangle(Point(0.2,0.2),Point(0.8,0.8))
#circle1 = Rectangle(Point(xinf,yinf),Point(0.15,0.15))
#circle2 = Rectangle(Point(0.85,0.15),Point(1.0,0))
#circle3 = Rectangle(Point(0.85,0.85),Point(xsup,ysup))
#circle4 = Rectangle(Point(0.15,0.85),Point(0,1.0))
#domain=rect-(circle+circle1+circle2+circle3+circle4)
#res = 80  # Resolution of mesh
#mesh = generate_mesh(domain, res)
#
#plot(mesh)
#plt.show()
#
#class Obstacle(SubDomain):
#    def inside(self, x, on_boundary):
#        return (on_boundary and between(x[0], (0.19,0.81)) and between(x[1],(0.19,0.81)))
#
#
#class Obstacle1(SubDomain):
#    def inside(self, x, on_boundary):
#        return (on_boundary and between(x[0], (-0.1,0.16)) and between(x[1], (-0.1,0.16)))
#
#
#class Obstacle2(SubDomain):
#    def inside(self, x, on_boundary):
#        return (on_boundary and between(x[0], (0.84,1.1)) and between(x[1], (-0.1,0.16)))
#
#
#class Obstacle3(SubDomain):
#    def inside(self, x, on_boundary):
#        return (on_boundary and between(x[0], (0.84,1.1)) and between(x[1], (0.84,1.1)))
#
#
#class Obstacle4(SubDomain):
#    def inside(self, x, on_boundary):
#        return (on_boundary and between(x[0], (-0.1,0.16)) and between(x[1], (0.84,1.1)))
#
#
#V=VectorFunctionSpace(mesh, 'P', 2, constrained_domain=PeriodicBoundary())
#nb_noeuds_ess=V.dim()
#
#print ('uuuuuuu',nb_noeuds_ess)
#
#obstacle = Obstacle()
#obstacle1 = Obstacle1()
#obstacle2 = Obstacle2()
#obstacle3 = Obstacle3()
#obstacle4 = Obstacle4()
#
#boundaries = FacetFunction("size_t", mesh)
#
#boundaries.set_all(5)
#obstacle.mark(boundaries, 0)
#obstacle1.mark(boundaries, 1)
#obstacle2.mark(boundaries, 2)
#obstacle3.mark(boundaries, 3)
#obstacle4.mark(boundaries, 4)
#
#    
#ds = Measure("ds")[boundaries]
#normale = FacetNormal(mesh) ## normale associe a chaque cellule
#
#num0=0
#num1=1
#num2=2
#num3=3
#num4=4
#
#
##################################################
######Sauvegarder les phi apres interpolation#####
##################################################  
#w=np.zeros((nb_noeuds_ess,nb_modes))
#wi_fixe=Function(V_fixe)
#
#repertoire_final_1="Homog_base_maill_init_multi/"
#
#File("Homog_base_maill_init_multi/mesh_complexe.xml.gz") << mesh
#
#ufile,USAVE=creation_fichier_pourecriture_champ_hdf5(repertoire_final_1,mesh)
#file_rayon_ecriture = open("%s/rayon_ecriture.txt" %(repertoire_final_1), "w")
#kfic=0
#for n in range(nb_modes):
#    xc=n+1
#    wi_fixe.vector().set_local(G[:,n])
#    wi_fixe.set_allow_extrapolation(True)
#    wi_ess=interpolate(wi_fixe,V)
##    plot(wi_ess)
##    plt.show()
#    ecriture_champ_hdf5(ufile,USAVE,wi_ess,kfic,file_rayon_ecriture,xc)
#    kfic=kfic+1
#    w[:,n]=wi_ess.vector().array()
#    
#
#
######Model reduit pour un rayon précis#####
#mat_mass=np.zeros((nb_modes,nb_modes))
#vec_sf=np.zeros((nb_modes))
#phiu_i = Function(V)                                                                    
#phiu_j = Function(V)                                                                    
#    
#
#
#for iproj in range(nb_modes):
#    phiu_i.vector().set_local(w[:,iproj])
#	
#        
#    for jrec in range(nb_modes):
#          phiu_j.vector().set_local(w[:,jrec])
#              
#              
#          ### matrice de masse
#          mat_mass[jrec,iproj]=assemble(tr(dot((grad(phiu_j)).T, grad(phiu_i)))*dx)
#          
#
#for iproj in range(nb_modes):
#    phiu_i.vector().set_local(w[:,iproj])
#    ### terme sf
#    vec_sf[iproj]=assemble((dot(normale,phiu_i)*ds(num0))+(dot(normale,phiu_i)*ds(num1))+(dot(normale,phiu_i)*ds(num2))+(dot(normale,phiu_i)*ds(num3))+(dot(normale,phiu_i)*ds(num4)))
#
#sol=np.linalg.solve(mat_mass.T,-vec_sf)
#
#print sol
#
#repertoire_final_2="Solution_homog_red_maill_init_multi/"
#
#ufile,USAVE=creation_fichier_pourecriture_champ_hdf5(repertoire_final_2,mesh)
#file_rayon_ecriture = open("%s/rayon_ecriture.txt" %(repertoire_final_2), "w")
#kfic=0
#
#
#urec=np.dot(w,sol)
#urec_red=Function(V)
#urec_red.vector().set_local(urec)
#
#ecriture_champ_hdf5(ufile,USAVE,urec_red,kfic,file_rayon_ecriture,1)
#
#H=assemble(grad(urec_red)[0,0]*dx)
#M=assemble(grad(urec_red)[1,1]*dx)
#X=assemble(grad(urec_red)[0,1]*dx)
#C=assemble(grad(urec_red)[1,0]*dx)
#print H,M,X,C
#ep=(1-0.3*0.3-4*0.15*0.15)
#D=((H/ep)+1)*ep
#print ("hi",r,ep,D)
#
#
#
#plot(urec_red)
#plt.show()
