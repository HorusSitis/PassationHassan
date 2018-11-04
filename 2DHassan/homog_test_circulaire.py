# -*- coding: utf-8 -*-
"""
Created on Fri Jul 13 17:57:52 2018

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

r=0.28



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
#circle = Circle(Point(c_x,c_y),r)
#domain=rect-circle
#res = 50  # Resolution of mesh
#mesh = generate_mesh(domain, res)
#
##### Premier raffinement dans l'axe horizontale du cylindre 
#meshB= raffinemment_maillage(r,mesh)
#mesh=meshB

mesh = Mesh("Homog_base_maill_init_circulaire/mesh_circulaire.xml.gz")

    
    
#class Obstacle(SubDomain):
#    def inside(self, x, on_boundary):
#        return (on_boundary and between((x[0]-c_x), (-r-tol, r+tol)) and between((x[1]-c_x), (-r-tol, r+tol)))
class Obstacle(SubDomain):
    def inside(self, x, on_boundary):
        return (on_boundary and between((x[0]-c_x), (-r-tol, r+tol)) and between((x[1]-c_x), (-r-tol, r+tol)))


V=VectorFunctionSpace(mesh, 'P', 2, constrained_domain=PeriodicBoundary())

ncel=V.dim()
obstacle = Obstacle()
    
P_impose=Constant((0., 0.))
pref = DirichletBC(V, P_impose, "x[0] < DOLFIN_EPS && x[1] < DOLFIN_EPS", "pointwise")
    
  

boundaries = FacetFunction("size_t", mesh)

boundaries.set_all(0)
obstacle.mark(boundaries, 5)

ds = Measure("ds")[boundaries]
normale = FacetNormal(mesh) ## normale associe a chaque cellule
    
num_front_cercle=5 ## numero de la frontiere du cylindre
    
    
u = TrialFunction(V)
v = TestFunction(V)
    
a=tr(dot((grad(u)).T, grad(v)))*dx
L=-dot(normale,v)*ds(num_front_cercle)
    
  
u=Function(V)
solve(a==L,u,pref)

plot(u)
plt.show()

H=assemble(grad(u)[0,0]*dx)
M=assemble(grad(u)[1,1]*dx)
X=assemble(grad(u)[0,1]*dx)
C=assemble(grad(u)[1,0]*dx)
print H,M,X,C
ep=(1-pi*(r**2))
D=((H/ep)+1)*ep
print ("hi",r,ep,D)

################################################################
###############test A[:,0] par psi*phi dans V_fixe##############
################################################################
#test=np.zeros((6))
#u.set_allow_extrapolation(True)
#u1=interpolate(u,V_fixe)
#repertoire_read="Homog_base_maill_init/"
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
#repertoire_read="Homog_base_maill_init_circulaire/"
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
repertoire_read="Solution_homog_red_maill_init_circulaire/"
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