# -*- coding: utf-8 -*-
"""
Created on Mon Jul  9 11:52:27 2018

@author: ghraieb
"""

from fenics import *
from dolfin import *
from mshr import *
import matplotlib.pyplot as plt
import numpy as np
from math import *
from Lecture_ecriture_homog import *
import sys
import os


xinf=0.0
yinf=0.0
xsup=1.0
ysup=1.0
c_x=0.5
c_y=0.5



def raffinemment_maillage(r,mesh):
    markers = CellFunction("bool", mesh)
    markers.set_all(False)
    for c in cells(mesh):
	# Mark cells with facet midpoints near y == 1.0
	for f in facets(c):
	    if (sqrt((f.midpoint()[0]-c_x)**2+(f.midpoint()[1]-c_y)**2)<=1.4*r):
		markers[c] = True
    #new_mesh=refine(mesh, markers, redistribute=False)  
    new_mesh=refine(mesh, markers, redistribute=True)  
    
    return new_mesh

R=np.array([0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.49,0.495])

domaine_fixe=Rectangle(Point(xinf,yinf),Point(xsup,ysup))
mesh_fixe=generate_mesh(domaine_fixe,50)
V_fixe=VectorFunctionSpace(mesh_fixe, "P", 2)


repertoire_final="Solutions_homog_fixe/"

ufile,USAVE=creation_fichier_pourecriture_champ_hdf5(repertoire_final,mesh_fixe)
file_rayon_ecriture = open("%s/rayon_ecriture.txt" %(repertoire_final), "w")
kfic=0

for i in(R):
    r=i

    rect= Rectangle(Point(xinf,yinf),Point(xsup,ysup))
    circle = Circle(Point(c_x,c_y),r)
    domain=rect-circle
    res = 25  # Resolution of mesh
    mesh = generate_mesh(domain, res)
#    plot(mesh)
#    plt.show()
    
    


    #### Premier raffinement dans l axe horizontale du cylindre 
    meshB= raffinemment_maillage(r,mesh)
    mesh=meshB
    print r,mesh
#    plot(mesh)
#    plt.show()
    
    
    tol=1e-14
    
    
    
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


            
         

    class Obstacle(SubDomain):
        def inside(self, x, on_boundary):
            return (on_boundary and between((x[0]-c_x), (-2*r-tol, 2*r+tol)) and between((x[1]-c_x), (-2*r-tol, 2*r+tol)))


    V=VectorFunctionSpace(mesh, 'P', 2, constrained_domain=PeriodicBoundary())

    
    obstacle = Obstacle()
    
    P_impose=Constant((0., 0.))
    pref = DirichletBC(V, P_impose, "x[0] < DOLFIN_EPS && x[1] < DOLFIN_EPS", "pointwise")
    
  

    boundaries = FacetFunction("size_t", mesh)

    boundaries.set_all(0)
    obstacle.mark(boundaries, 5)
    
    ds = Measure("ds")[boundaries]
    normale = FacetNormal(mesh) ## normale associe a chaque cellule
    
    num_front_cercle=5 ## numero de la frontiere du cylindre
    
    nb_noeuds=V.dim()
    print nb_noeuds
    
    u = TrialFunction(V)
    v = TestFunction(V)
    
    a=tr(dot((grad(u)).T, grad(v)))*dx
    L=-dot(normale,v)*ds(num_front_cercle)
    
#    A=assemble(a)
#    B=assemble(L)
    
    u=Function(V)
#    bc.apply(A, B)
#    solve(A, u.vector(), B)
    solve(a==L,u,pref)
    
#    print u(0,0),u(1,1),u(0,1),u(1,0),u(0,0.5)-u(1,0.5),u(0.5,0)-u(0.5,1),u(0,0.6),u(1,0.6),u(0.4,0),u(0.4,1)

    H=assemble(grad(u)[0,0]*dx)
    M=assemble(grad(u)[1,1]*dx)
    A=assemble(grad(u)[0,1]*dx)
    C=assemble(grad(u)[1,0]*dx)
    print H,M,A,C
    ep=(1-pi*(r**2))
    D=((H/ep)+1)*ep
    print ("hi",r,ep,D)
    
    plot(u)
    plt.show()
    

#    parameters["allow_extrapolation"] = True    
    
    u.set_allow_extrapolation(True)
    u_fixe=interpolate(u,V_fixe)
    plot(u_fixe)
    plt.show()
    
    ecriture_champ_hdf5(ufile,USAVE,u_fixe,kfic,file_rayon_ecriture,r)
    kfic=kfic+1
    
    
    
    
    
    

    
        
    
    
    

            
