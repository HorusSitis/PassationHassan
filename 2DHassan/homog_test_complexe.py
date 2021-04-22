# -*- coding: utf-8 -*-
"""
Created on Tue Jul 17 16:06:27 2018

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

#R=np.arange(0.01,0.26,0.005)
R=[0.01,0.05,0.1,0.15,0.2,0.22,0.24,0.25,0.255]

r=0.232

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

domaine_fixe=Rectangle(Point(xinf,yinf),Point(xsup,ysup))
mesh_fixe=generate_mesh(domaine_fixe,100)
V_fixe=VectorFunctionSpace(mesh_fixe, "P", 2, constrained_domain=PeriodicBoundary())

#rect= Rectangle(Point(xinf,yinf),Point(xsup,ysup))
#circle = Circle(Point(c_x,c_y),0.45)
#circle1 = Circle(Point(c1_x,c1_y),r)
#circle2 = Circle(Point(c2_x,c2_y),r)
#circle3 = Circle(Point(c3_x,c3_y),r)
#circle4 = Circle(Point(c4_x,c4_y),r)
#domain=rect-(circle+circle1+circle2+circle3+circle4)
#res = 50  # Resolution of mesh
#mesh = generate_mesh(domain, res)
##### Premier raffinement dans l axe horizontale du cyٍlindre 
#meshB= raffinemment_maillage(r,mesh)
#mesh=meshB


mesh = Mesh("Homog_base_maill_init_complexe/mesh_complexe.xml.gz")



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
ncel=V.dim()
print ncel

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

u = TrialFunction(V)
v = TestFunction(V)

a=tr(dot((grad(u)).T, grad(v)))*dx
L=-dot(normale,v)*ds(num0)-dot(normale,v)*ds(num1)-dot(normale,v)*ds(num2)-dot(normale,v)*ds(num3)-dot(normale,v)*ds(num4)

u=Function(V)

solve(a==L,u)

H=assemble(grad(u)[0,0]*dx)
M=assemble(grad(u)[1,1]*dx)
A=assemble(grad(u)[0,1]*dx)
C=assemble(grad(u)[1,0]*dx)
print H,M,A,C
ep=(1-pi*((r**2)+0.45**2))
D=((H/ep)+1)*ep
print ("hi",r,ep,D)
plot(u)
plt.show()
################################################################
###############test A[:,0] par psi*phi dans V_fixe##############
################################################################
#test=np.zeros((6))
#u.set_allow_extrapolation(True)
#u1=interpolate(u,V_fixe)
#repertoire_read="Homog_base_maill_init_complexe/"
#for n in range(6):
#    w = Function(V)
#    xc=n+1
#    # Update current raduis
#    lecture_fichiers_restart_hdf5(w,repertoire_read,xc,mesh) 
#    w.set_allow_extrapolation(True)
#    u2=interpolate(w,V_fixe)
#    test[n]=assemble(dot(u1,u2)*dx)

##################################################################
#####################test A[:,0] par psi*phi dans V ##############
##################################################################

#u1=u.vector().array()
#test=np.zeros((6))
#repertoire_read="Homog_base_maill_init_complexe/"
#for n in range(6):
#    w = Function(V)
#    xc=n+1
#    # Update current raduis
#    lecture_fichiers_restart_hdf5(w,repertoire_read,xc,mesh) 
#    test[n]=assemble(dot(u,w)*dx)




##############################################
##################calcul erreur ##############
##############################################

ured=Function(V)
repertoire_read="Solution_homog_red_maill_init_complexe/"
lecture_solution_restart_hdf5(ured,repertoire_read,1,mesh)
plot(ured)
plt.show()

num=Function(V)
den=Function(V)

u1=u.vector().array()
u2=ured.vector().array()
u3=u1-u2
num.vector().set_local(u3)
den.vector().set_local(u1)

erreur=(sqrt(assemble(dot(num,num)*dx)))/(sqrt(assemble(dot(den,den)*dx)))

print erreur