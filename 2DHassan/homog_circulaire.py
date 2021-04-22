# -*- coding: utf-8 -*-
"""
Created on Wed Aug 29 15:45:25 2018

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



def raffinemment_maillage(r,mesh):
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

R=np.arange(0.1,0.5,0.05)
#R=np.zeros((10))
#R[0:8]=np.arange(0.1,0.5,0.05)
#R[8]=0.49
#R[9]=0.495

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


#determiner le domaine fixe pour interpoler la solution
domaine_fixe=Rectangle(Point(xinf,yinf),Point(xsup,ysup))
mesh_fixe=generate_mesh(domaine_fixe,80)
V_fixe=VectorFunctionSpace(mesh_fixe, "P", 2, constrained_domain=PeriodicBoundary())

plot(mesh_fixe)
plt.show()

repertoire_final="Solutions_homog_interp_circulaire/"

File("Solutions_homog_interp_circulaire/mesh_circulaire.xml.gz") << mesh_fixe


ufile,USAVE=creation_fichier_pourecriture_champ_hdf5(repertoire_final,mesh_fixe)
file_rayon_ecriture = open("%s/rayon_ecriture.txt" %(repertoire_final), "w")
kfic=0


for i in R:
    r=i

    rect= Rectangle(Point(xinf,yinf),Point(xsup,ysup))
    circle = Circle(Point(c_x,c_y),r)
    domain=rect-circle
    res = 40  # Resolution of mesh
    mesh = generate_mesh(domain, res)

    #### Premier raffinement dans l'axe horizontale du cylindre 
    meshB= raffinemment_maillage(r,mesh)
    mesh=meshB
    plot(mesh)
    plt.show()
    
    
#    class Obstacle(SubDomain):
#        def inside(self, x, on_boundary):
#            return (on_boundary and between((x[0]-c_x), (-r-tol, r+tol)) and between((x[1]-c_x), (-r-tol, r+tol)))
    class Obstacle(SubDomain):
        def inside(self, x, on_boundary):
            return (on_boundary and between((x[0]-c_x), (-r-tol, r+tol)) and between((x[1]-c_x), (-r-tol, r+tol)))


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
#    print nb_noeuds
    
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
    
#    print u(0,0.5),u(c_x+r,c_y),u(c_x,c_y-r),u(c_x,c_y+r)#u(0.1,0.1),u(0.1,0.9),u(0.9,0.1),u(0.9,0.9),u(0,0.5),u(1,0.5),u(0.5,0),u(0.5,1),u(0.2,0.2),u(0,0.6),u(1,0.6),u(0.4,0),u(0.4,1)

    H=assemble(grad(u)[0,0]*dx)
    M=assemble(grad(u)[1,1]*dx)
    A=assemble(grad(u)[0,1]*dx)
    C=assemble(grad(u)[1,0]*dx)
    print H,M,A,C
    ep=(1-pi*(r**2))
    D=((H/ep)+1)*ep
    print ("hi",r,ep,D)
    
#    if r==R[0] or r==R[3] or r==R[7]:
    plot(u)
    plt.show()
    
    u.set_allow_extrapolation(True)
    u_fixe1=interpolate(u,V_fixe)
#    u_fixe1 = project(u, V_fixe)
        
    plot(u_fixe1)
    plt.show()

    
#    val1=u_fixe1.vector().array()
#    
#    fi=0.0
#    fe=1.0
#    fcar=Expression(("sqrt(pow((x[0]-0.5),2)+pow((x[1]-0.5),2)) - r < 0.0 ? fi : fe", "sqrt(pow((x[0]-0.5),2)+pow((x[1]-0.5),2)) - r < 0.0 ? fi : fe"),degree=2, fi=fi, fe=fe, r=r)
#    
#    u_car=interpolate(fcar,V_fixe)
#    val2=u_car.vector().array()
#    
#    d_dim=len(val1)
#    val3=np.zeros((d_dim))
#    for n in range(d_dim):
#        val3[n]=val1[n]*val2[n]
#    
#    u_fixe=Function(V_fixe)
#    u_fixe.vector().set_local(val3)
#    
#
##    print u_fixe(0,0.5),u_fixe(c_x+r,c_y),u_fixe(c_x,c_y-r),u_fixe(c_x,c_y+r),u_fixe(0.5,0.5),u_fixe(0.45,0.45)
#
#    ############################################
#    #########test coef homog####################
#    ############################################
#
#    class Omega_fluide(SubDomain):
#        def inside(self, x, on_boundary):
#            return (((sqrt((x[0] - c_x)**2 + (x[1] - c_y)**2)-r) >= 0.0))
#
#    subdomain_fluide = Omega_fluide()
#    subdomains = CellFunction('size_t', mesh_fixe)
#    default_value=0
#    subdomains.set_all(default_value)    
#    subdomain_fluide.mark(subdomains, 1)
#    dx1 = Measure("dx")[subdomains]
#    num_surf=1    
#    H=assemble(grad(u_fixe)[0,0]*dx1(num_surf))
#    M=assemble(grad(u_fixe)[1,1]*dx1(num_surf))
#    A=assemble(grad(u_fixe)[0,1]*dx1(num_surf))
#    C=assemble(grad(u_fixe)[1,0]*dx1(num_surf))
#    print H,M,A,C
#    ep=(1-pi*(r**2))
#    D=((H/ep)+1)*ep
#    print ("hi",r,ep,D)
    
#    print u_fixe(0,0.5),u_fixe(1,0.5),u_fixe(0.5,0),u_fixe(0.5,1),u_fixe(0,0.4),u_fixe(1,0.4),u_fixe(0.5,0.5),u_fixe(1-r-0.01,0),u_fixe(1-r-0.01,1)

    
#    plot(u_fixe)
#    plt.show()

#    if r==R[0] or r==R[3] or r==R[7]:
#        plot(u_fixe)
#        plt.show()
    
    ecriture_champ_hdf5(ufile,USAVE,u_fixe1,kfic,file_rayon_ecriture,r)
    kfic=kfic+1