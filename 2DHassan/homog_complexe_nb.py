# -*- coding: utf-8 -*-
"""
Created on Tue Aug 28 15:29:53 2018

@author: ghraieb
"""

from fenics import *
from dolfin import *
from mshr import *
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

#R=np.arange(0.01,0.27,0.02)
R=np.arange(0.01,0.26,0.01)

def raffinemment_maillage(r,mesh):
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
mesh_fixe=generate_mesh(domaine_fixe,100)
V_fixe=VectorFunctionSpace(mesh_fixe, "P", 2, constrained_domain=PeriodicBoundary())

plot(mesh_fixe)
plt.show()

repertoire_final="Solutions_homog_interp_complexe_nb/"

File("Solutions_homog_interp_complexe_nb/mesh_complexe_nb.xml.gz") << mesh_fixe

ufile,USAVE=creation_fichier_pourecriture_champ_hdf5(repertoire_final,mesh_fixe)
file_rayon_ecriture = open("%s/rayon_ecriture.txt" %(repertoire_final), "w")
kfic=0


for i in(R):
    r=i
    rect= Rectangle(Point(xinf,yinf),Point(xsup,ysup))
    circle = Circle(Point(c_x,c_y),0.45)
    circle1 = Circle(Point(c1_x,c1_y),r)
    circle2 = Circle(Point(c2_x,c2_y),r)
    circle3 = Circle(Point(c3_x,c3_y),r)
    circle4 = Circle(Point(c4_x,c4_y),r)
    domain=rect-(circle+circle1+circle2+circle3+circle4)
    res = 50  # Resolution of mesh
    mesh = generate_mesh(domain, res)

    #### Premier raffinement dans l axe horizontale du cyٍlindre 
    meshB= raffinemment_maillage(r,mesh)
    mesh=meshB
#    plot(mesh)
#    plt.show()


    class Obstacle(SubDomain):
        def inside(self, x, on_boundary):
            return (on_boundary and between((x[0]-c_x), (-0.45-tol, 0.45+tol)) and between((x[1]-c_x), (-0.45-tol, 0.45+tol)))


    class Obstacle1(SubDomain):
        def inside(self, x, on_boundary):
            return (on_boundary and between((x[0]-c1_x), (-r-tol, r+tol)) and between((x[1]-c1_y), (-r-tol, r+tol)))


    class Obstacle2(SubDomain):
        def inside(self, x, on_boundary):
            return (on_boundary and between((x[0]-c2_x), (-r-tol, r+tol)) and between((x[1]-c2_y), (-r-tol, r+tol)))


    class Obstacle3(SubDomain):
        def inside(self, x, on_boundary):
            return (on_boundary and between((x[0]-c3_x), (-r-tol, r+tol)) and between((x[1]-c3_y), (-r-tol, r+tol)))


    class Obstacle4(SubDomain):
        def inside(self, x, on_boundary):
            return (on_boundary and between((x[0]-c4_x), (-r-tol, r+tol)) and between((x[1]-c4_y), (-r-tol, r+tol)))


    V=VectorFunctionSpace(mesh, 'P', 2, constrained_domain=PeriodicBoundary())

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

    nb_noeuds=V.dim()
    
    u = TrialFunction(V)
    v = TestFunction(V)
    
    a=tr(dot((grad(u)).T, grad(v)))*dx
    L=-dot(normale,v)*ds(num0)-dot(normale,v)*ds(num1)-dot(normale,v)*ds(num2)-dot(normale,v)*ds(num3)-dot(normale,v)*ds(num4)

    u=Function(V)

    solve(a==L,u)

#    print u(0.1,0),u(0.1,1),u(0,0.1),u(1,0.1),u(0.5,0.5)#u(0,0.5),u(1,0.5),u(0.5,0),u(0.5,1),u(0,0.4),u(1,0.4),u(0.5+0.45,0.5),u(1-r-0.01,0),u(1-r-0.01,1)
#    if r==R[4]:
#        plot(u)
#        plt.show()
    
    H=assemble(grad(u)[0,0]*dx)
    M=assemble(grad(u)[1,1]*dx)
    A=assemble(grad(u)[0,1]*dx)
    C=assemble(grad(u)[1,0]*dx)
    print H,M,A,C
    ep=(1-pi*((r**2)+0.45**2))
    D=((H/ep)+1)*ep
    print ("hi",r,ep,D)

#    plot(u)
#    plt.show()

    u.set_allow_extrapolation(True)
    u_fixe1=interpolate(u,V_fixe)

    
    val1=u_fixe1.vector().array()
    
    fi=0.0
    fe=1.0
    fcar=Expression(("((sqrt(pow((x[0]-0.5),2)+pow((x[1]-0.5),2)) - 0.45 < 0.0) or (sqrt(pow((x[0]-c4_x),2)+pow((x[1]-c4_y),2)) - r < 0.0) or (sqrt(pow((x[0]-c1_x),2)+pow((x[1]-c1_y),2)) - r < 0.0) or (sqrt(pow((x[0]-c2_x),2)+pow((x[1]-c2_y),2)) - r < 0.0) or (sqrt(pow((x[0]-c3_x),2)+pow((x[1]-c3_y),2)) - r < 0.0)) ? fi : fe", "((sqrt(pow((x[0]-0.5),2)+pow((x[1]-0.5),2)) - 0.45 < 0.0) or (sqrt(pow((x[0]-c1_x),2)+pow((x[1]-c1_y),2)) - r < 0.0) or (sqrt(pow((x[0]-c2_x),2)+pow((x[1]-c2_y),2)) - r < 0.0) or (sqrt(pow((x[0]-c3_x),2)+pow((x[1]-c3_y),2)) - r < 0.0) or (sqrt(pow((x[0]-c4_x),2)+pow((x[1]-c4_y),2)) - r < 0.0)) ? fi : fe"),degree=2, fi=fi, fe=fe, r=r,c1_x=c1_x,c1_y=c1_y,c2_x=c2_x,c2_y=c2_y,c3_x=c3_x,c3_y=c3_y,c4_x=c4_x,c4_y=c4_y)
    
    u_car=interpolate(fcar,V_fixe)
    val2=u_car.vector().array()
    
    d_dim=len(val1)
    val3=np.zeros((d_dim))
    for n in range(d_dim):
        val3[n]=val1[n]*val2[n]
    
    u_fixe=Function(V_fixe)
    u_fixe.vector().set_local(val3)

#    plot(u_fixe)
#    plt.show()
    
    ############################################
    #########test coef homog####################
    ############################################

    class Omega_fluide(SubDomain):
        def inside(self, x, on_boundary):
            return (((sqrt((x[0] - c_x)**2 + (x[1] - c_y)**2)-0.45) >= 0.0) and ((sqrt((x[0] - c1_x)**2 + (x[1] - c1_y)**2)-r) >= 0.0) and ((sqrt((x[0] - c2_x)**2 + (x[1] - c2_y)**2)-r) >= 0.0) and ((sqrt((x[0] - c3_x)**2 + (x[1] - c3_y)**2)-r) >= 0.0) and ((sqrt((x[0] - c4_x)**2 + (x[1] - c4_y)**2)-r) >= 0.0))

    subdomain_fluide = Omega_fluide()
    subdomains = CellFunction('size_t', mesh_fixe)
    default_value=0
    subdomains.set_all(default_value)    
    subdomain_fluide.mark(subdomains, 1)
    dx1 = Measure("dx")[subdomains]
    num_surf=1    
    H=assemble(grad(u_fixe)[0,0]*dx1(num_surf))
    M=assemble(grad(u_fixe)[1,1]*dx1(num_surf))
    A=assemble(grad(u_fixe)[0,1]*dx1(num_surf))
    C=assemble(grad(u_fixe)[1,0]*dx1(num_surf))
    print H,M,A,C
    ep=(1-pi*((r**2)+0.45**2))
    D=((H/ep)+1)*ep
    print ("hi",r,ep,D)

#    print u_fixe(0.005,0),u_fixe(0.005,1),u_fixe(0,0.005),u_fixe(1,0.005),u_fixe(0.5,0.5)#u(0,0.5),u(1,0.5),u(0.5,0),u(0.5,1),u(0,0.4),u(1,0.4),u(0.5+0.45,0.5),u(1-r-0.01,0),u(1-r-0.01,1)
    
#    print u_fixe(0,0.5),u_fixe(1,0.5),u_fixe(0.5,0),u_fixe(0.5,1),u_fixe(0,0.4),u_fixe(1,0.4),u_fixe(0.5+0.45,0.5),u_fixe(1-r-0.01,0),u_fixe(1-r-0.01,1),u_fixe(0.5,0.5)

    
#    plot(u_fixe)
#    plt.show()
    
    ecriture_champ_hdf5(ufile,USAVE,u_fixe1,kfic,file_rayon_ecriture,r)
    kfic=kfic+1

#    if r==R[4]:
#        break