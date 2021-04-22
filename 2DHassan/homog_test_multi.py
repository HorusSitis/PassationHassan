# -*- coding: utf-8 -*-
"""
Created on Thu Sep  6 08:28:46 2018

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


###############################################
###################circulaire##################
###############################################
###################test réduit#################
###############################################
#r=0.423
#
#
#mesh = Mesh("Homog_base_maill_init_multi/mesh_circulaire.xml.gz")
#
#class Obstacle(SubDomain):
#    def inside(self, x, on_boundary):
#        return (on_boundary and between((x[0]-c_x), (-r-tol, r+tol)) and between((x[1]-c_x), (-r-tol, r+tol)))
#
#
#V=VectorFunctionSpace(mesh, 'P', 2, constrained_domain=PeriodicBoundary())
#
#ncel=V.dim()
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
#    
#u = TrialFunction(V)
#v = TestFunction(V)
#    
#a=tr(dot((grad(u)).T, grad(v)))*dx
#L=-dot(normale,v)*ds(num_front_cercle)
#    
#  
#u=Function(V)
#solve(a==L,u,pref)
#
#plot(u)
#plt.show()
#
#H=assemble(grad(u)[0,0]*dx)
#M=assemble(grad(u)[1,1]*dx)
#X=assemble(grad(u)[0,1]*dx)
#C=assemble(grad(u)[1,0]*dx)
#print H,M,X,C
#ep=(1-pi*(r**2))
#D=((H/ep)+1)*ep
#print ("hi",r,ep,D)
#
###############################################
###################calcul erreur ##############
###############################################
#
#ured=Function(V)
#repertoire_read="Solution_homog_red_maill_init_multi/"
#lecture_solution_restart_hdf5(ured,repertoire_read,1,mesh)
#plot(ured)
#plt.show()
#
#num=Function(V)
#den=Function(V)
#
#u1=u.vector().array()
#u2=ured.vector().array()
#u3=u1-u2
#num.vector().set_local(u3)
#den.vector().set_local(u1)
#
#erreur=(sqrt(assemble(dot(num,num)*dx)))/(sqrt(assemble(dot(den,den)*dx)))
#
#print erreur





###############################################
###################complexe_1##################
###############################################
###################test réduit#################
###############################################
#r=0.087
#mesh = Mesh("Homog_base_maill_init_multi/mesh_complexe_1.xml.gz")
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
#ncel=V.dim()
#
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
#ds = Measure("ds")[boundaries]
#normale = FacetNormal(mesh) ## normale associe a chaque cellule
#
#
#num1=1
#num2=2
#num3=3
#num4=4
#
#u = TrialFunction(V)
#v = TestFunction(V)
#
#a=tr(dot((grad(u)).T, grad(v)))*dx
#L=-dot(normale,v)*ds(num1)-dot(normale,v)*ds(num2)-dot(normale,v)*ds(num3)-dot(normale,v)*ds(num4)
#
#u=Function(V)
#
#solve(a==L,u)
#
#H=assemble(grad(u)[0,0]*dx)
#M=assemble(grad(u)[1,1]*dx)
#A=assemble(grad(u)[0,1]*dx)
#C=assemble(grad(u)[1,0]*dx)
#print H,M,A,C
#ep=(1-pi*(r**2))
#D=((H/ep)+1)*ep
#print ("hi",r,ep,D)
#
#plot(u)
#plt.show()
#
###############################################
###################calcul erreur ##############
###############################################
#
#ured=Function(V)
#repertoire_read="Solution_homog_red_maill_init_multi/"
#lecture_solution_restart_hdf5(ured,repertoire_read,1,mesh)
#plot(ured)
#plt.show()
#
#num=Function(V)
#den=Function(V)
#err=Function(V)
#
#u1=u.vector().array()
#u2=ured.vector().array()
#u3=u1-u2
#u4=(u3/u1)*100
##for i in range(3*ncel):
##    print u1[i],u2[i],u4[i]
#err.vector().set_local(u4)
#num.vector().set_local(u3)
#den.vector().set_local(u1)
#plot(err)
#plt.show()
#erreur=(sqrt(assemble(dot(num,num)*dx)))/(sqrt(assemble(dot(den,den)*dx)))
#
#print erreur





###############################################
####################complexe###################
###############################################
###################test réduit#################
###############################################
#r=0.248
#mesh = Mesh("Homog_base_maill_init_multi/mesh_complexe.xml.gz")
#
#
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
#ncel=V.dim()
#print ncel
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
#ds = Measure("ds")[boundaries]
#normale = FacetNormal(mesh) ## normale associe a chaque cellule
#
#num0=0
#num1=1
#num2=2
#num3=3
#num4=4
#
#u = TrialFunction(V)
#v = TestFunction(V)
#
#a=tr(dot((grad(u)).T, grad(v)))*dx
#L=-dot(normale,v)*ds(num0)-dot(normale,v)*ds(num1)-dot(normale,v)*ds(num2)-dot(normale,v)*ds(num3)-dot(normale,v)*ds(num4)
#
#u=Function(V)
#
#solve(a==L,u)
#
#H=assemble(grad(u)[0,0]*dx)
#M=assemble(grad(u)[1,1]*dx)
#A=assemble(grad(u)[0,1]*dx)
#C=assemble(grad(u)[1,0]*dx)
#print H,M,A,C
#ep=(1-pi*((r**2)+0.45**2))
#D=((H/ep)+1)*ep
#print ("hi",r,ep,D)
#plot(u)
#plt.show()
#
###############################################
###################calcul erreur ##############
###############################################
#
#ured=Function(V)
#repertoire_read="Solution_homog_red_maill_init_multi/"
#lecture_solution_restart_hdf5(ured,repertoire_read,1,mesh)
#plot(ured)
#plt.show()
#
#num=Function(V)
#den=Function(V)
#
#u1=u.vector().array()
#u2=ured.vector().array()
#u3=u1-u2
#num.vector().set_local(u3)
#den.vector().set_local(u1)
#
#erreur=(sqrt(assemble(dot(num,num)*dx)))/(sqrt(assemble(dot(den,den)*dx)))
#
#print erreur







###############################################
###########un quart et cercle central##########
###############################################
###################test réduit#################
###############################################
#r=0.2
#
#mesh = Mesh("Homog_base_maill_init_multi/mesh_quelconque.xml.gz")
#
#
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
#ncel=V.dim()
#print ncel
#
#obstacle = Obstacle()
#obstacle1 = Obstacle1()
#
#
#boundaries = FacetFunction("size_t", mesh)
#
#boundaries.set_all(5)
#obstacle.mark(boundaries, 0)
#obstacle1.mark(boundaries, 1)
#
#ds = Measure("ds")[boundaries]
#normale = FacetNormal(mesh) ## normale associe a chaque cellule
#
#num0=0
#num1=1
#
#u = TrialFunction(V)
#v = TestFunction(V)
#
#a=tr(dot((grad(u)).T, grad(v)))*dx
#L=-dot(normale,v)*ds(num0)-dot(normale,v)*ds(num1)
#
#u=Function(V)
#
#solve(a==L,u)
#
#print u(0.6,0),u(0.6,1),u(0,0.6),u(1,0.6),u(0,0.7),u(1,0.7)
#H=assemble(grad(u)[0,0]*dx)
#M=assemble(grad(u)[1,1]*dx)
#A=assemble(grad(u)[0,1]*dx)
#C=assemble(grad(u)[1,0]*dx)
#print H,M,A,C
#ep=(1-pi*(0.25*(r**2)+0.45**2))
#D=((H/ep)+1)*ep
#print ("hi",r,ep,D)
#plot(u)
#plt.show()
#
###############################################
###################calcul erreur ##############
###############################################
#
#ured=Function(V)
#repertoire_read="Solution_homog_red_maill_init_multi/"
#lecture_solution_restart_hdf5(ured,repertoire_read,1,mesh)
#plot(ured)
#plt.show()
#
#num=Function(V)
#den=Function(V)
#
#u1=u.vector().array()
#u2=ured.vector().array()
#u3=u1-u2
#num.vector().set_local(u3)
#den.vector().set_local(u1)
#
#erreur=(sqrt(assemble(dot(num,num)*dx)))/(sqrt(assemble(dot(den,den)*dx)))
#
#print erreur










###############################################
########### 3 quart et cercle central##########
###############################################
###################test réduit#################
###############################################
#r=0.2
#
#mesh = Mesh("Homog_base_maill_init_multi/mesh_quelconque.xml.gz")
#
#
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
#ncel=V.dim()
#print ncel
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
#ds = Measure("ds")[boundaries]
#normale = FacetNormal(mesh) ## normale associe a chaque cellule
#
#num0=0
#num1=1
#num2=2
#num3=3
#
#u = TrialFunction(V)
#v = TestFunction(V)
#
#a=tr(dot((grad(u)).T, grad(v)))*dx
#L=-dot(normale,v)*ds(num0)-dot(normale,v)*ds(num1)-dot(normale,v)*ds(num2)-dot(normale,v)*ds(num3)
#
#u=Function(V)
#
#solve(a==L,u)
#
#print u(0.6,0),u(0.6,1),u(0,0.6),u(1,0.6)
#H=assemble(grad(u)[0,0]*dx)
#M=assemble(grad(u)[1,1]*dx)
#A=assemble(grad(u)[0,1]*dx)
#C=assemble(grad(u)[1,0]*dx)
#print H,M,A,C
#ep=(1-pi*(0.75*(r**2)+0.45**2))
#D=((H/ep)+1)*ep
#print ("hi",r,ep,D)
#plot(u)
#plt.show()
#
###############################################
###################calcul erreur ##############
###############################################
#
#ured=Function(V)
#repertoire_read="Solution_homog_red_maill_init_multi/"
#lecture_solution_restart_hdf5(ured,repertoire_read,1,mesh)
#plot(ured)
#plt.show()
#
#num=Function(V)
#den=Function(V)
#
#u1=u.vector().array()
#u2=ured.vector().array()
#u3=u1-u2
#num.vector().set_local(u3)
#den.vector().set_local(u1)
#
#erreur=(sqrt(assemble(dot(num,num)*dx)))/(sqrt(assemble(dot(den,den)*dx)))
#
#print erreur










###############################################
#####2 quart sur même axe et cercle central####
###############################################
###################test réduit#################
###############################################
#r=0.2
#
#mesh = Mesh("Homog_base_maill_init_multi/mesh_quelconque.xml.gz")
#
#
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
#ncel=V.dim()
#print ncel
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
#ds = Measure("ds")[boundaries]
#normale = FacetNormal(mesh) ## normale associe a chaque cellule
#
#num0=0
#num1=1
#num2=2
#
#u = TrialFunction(V)
#v = TestFunction(V)
#
#a=tr(dot((grad(u)).T, grad(v)))*dx
#L=-dot(normale,v)*ds(num0)-dot(normale,v)*ds(num1)-dot(normale,v)*ds(num2)
#
#u=Function(V)
#
#solve(a==L,u)
#
#print u(0.6,0),u(0.6,1),u(0,0.6),u(1,0.6)
#H=assemble(grad(u)[0,0]*dx)
#M=assemble(grad(u)[1,1]*dx)
#A=assemble(grad(u)[0,1]*dx)
#C=assemble(grad(u)[1,0]*dx)
#print H,M,A,C
#ep=(1-pi*(0.5*(r**2)+0.45**2))
#D=((H/ep)+1)*ep
#print ("hi",r,ep,D)
#plot(u)
#plt.show()
#
###############################################
###################calcul erreur ##############
###############################################
#
#ured=Function(V)
#repertoire_read="Solution_homog_red_maill_init_multi/"
#lecture_solution_restart_hdf5(ured,repertoire_read,1,mesh)
#plot(ured)
#plt.show()
#
#num=Function(V)
#den=Function(V)
#
#u1=u.vector().array()
#u2=ured.vector().array()
#u3=u1-u2
#num.vector().set_local(u3)
#den.vector().set_local(u1)
#
#erreur=(sqrt(assemble(dot(num,num)*dx)))/(sqrt(assemble(dot(den,den)*dx)))
#
#print erreur










###############################################
#### 2 quart sur diagonal et cercle central####
###############################################
###################test réduit#################
###############################################
#r=0.2
#
#mesh = Mesh("Homog_base_maill_init_multi/mesh_quelconque.xml.gz")
#
#
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
#V=VectorFunctionSpace(mesh, 'P', 2, constrained_domain=PeriodicBoundary())
#ncel=V.dim()
#print ncel
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
#ds = Measure("ds")[boundaries]
#normale = FacetNormal(mesh) ## normale associe a chaque cellule
#
#num0=0
#num1=1
#num3=3
#
#u = TrialFunction(V)
#v = TestFunction(V)
#
#a=tr(dot((grad(u)).T, grad(v)))*dx
#L=-dot(normale,v)*ds(num0)-dot(normale,v)*ds(num1)-dot(normale,v)*ds(num3)
#
#u=Function(V)
#
#solve(a==L,u)
#
#print u(0.6,0),u(0.6,1),u(0,0.6),u(1,0.6)
#H=assemble(grad(u)[0,0]*dx)
#M=assemble(grad(u)[1,1]*dx)
#A=assemble(grad(u)[0,1]*dx)
#C=assemble(grad(u)[1,0]*dx)
#print H,M,A,C
#ep=(1-pi*(0.5*(r**2)+0.45**2))
#D=((H/ep)+1)*ep
#print ("hi",r,ep,D)
#plot(u)
#plt.show()
#
###############################################
###################calcul erreur ##############
###############################################
#
#ured=Function(V)
#repertoire_read="Solution_homog_red_maill_init_multi/"
#lecture_solution_restart_hdf5(ured,repertoire_read,1,mesh)
#plot(ured)
#plt.show()
#
#num=Function(V)
#den=Function(V)
#
#u1=u.vector().array()
#u2=ured.vector().array()
#u3=u1-u2
#num.vector().set_local(u3)
#den.vector().set_local(u1)
#
#erreur=(sqrt(assemble(dot(num,num)*dx)))/(sqrt(assemble(dot(den,den)*dx)))
#
#print erreur









###############################################
######4 differents quart et cercle central#####
###############################################
###################test réduit#################
###############################################
#r1=0.05
#r2=0.1
#r3=0.15
#r4=0.2
##r1=r2=r3=r4
#
#mesh = Mesh("Homog_base_maill_init_multi/mesh_complexe.xml.gz")
#
#
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
#ncel=V.dim()
#print ncel
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
#ds = Measure("ds")[boundaries]
#normale = FacetNormal(mesh) ## normale associe a chaque cellule
#
#num0=0
#num1=1
#num2=2
#num3=3
#num4=4
#
#u = TrialFunction(V)
#v = TestFunction(V)
#
#a=tr(dot((grad(u)).T, grad(v)))*dx
#L=-dot(normale,v)*ds(num0)-dot(normale,v)*ds(num1)-dot(normale,v)*ds(num2)-dot(normale,v)*ds(num3)-dot(normale,v)*ds(num4)
#
#u=Function(V)
#
#solve(a==L,u)
#
#H=assemble(grad(u)[0,0]*dx)
#M=assemble(grad(u)[1,1]*dx)
#A=assemble(grad(u)[0,1]*dx)
#C=assemble(grad(u)[1,0]*dx)
#print H,M,A,C
#ep=(1-pi*(0.25*(r1**2)+0.25*(r2**2)+0.25*(r3**2)+0.25*(r4**2)+0.45**2))
#D=((H/ep)+1)*ep
#print ("hi",ep,D)
#plot(u)
#plt.show()
#
###############################################
###################calcul erreur ##############
###############################################
#
#ured=Function(V)
#repertoire_read="Solution_homog_red_maill_init_multi/"
#lecture_solution_restart_hdf5(ured,repertoire_read,1,mesh)
#plot(ured)
#plt.show()
#
#num=Function(V)
#den=Function(V)
#
#u1=u.vector().array()
#u2=ured.vector().array()
#u3=u1-u2
#num.vector().set_local(u3)
#den.vector().set_local(u1)
#
#erreur=(sqrt(assemble(dot(num,num)*dx)))/(sqrt(assemble(dot(den,den)*dx)))
#
#print erreur











###############################################
############ 4 différentes quarts##############
###############################################
###################test réduit#################
###############################################
#r1=0.1
#r2=0.2
#r3=0.3
#r4=0.4
#mesh = Mesh("Homog_base_maill_init_multi/mesh_complexe.xml.gz")
#
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
#ncel=V.dim()
#print ncel
#
#obstacle1 = Obstacle1()
#obstacle2 = Obstacle2()
#obstacle3 = Obstacle3()
#obstacle4 = Obstacle4()
#
#boundaries = FacetFunction("size_t", mesh)
#
#boundaries.set_all(5)
#obstacle1.mark(boundaries, 1)
#obstacle2.mark(boundaries, 2)
#obstacle3.mark(boundaries, 3)
#obstacle4.mark(boundaries, 4)
#
#ds = Measure("ds")[boundaries]
#normale = FacetNormal(mesh) ## normale associe a chaque cellule
#
#num1=1
#num2=2
#num3=3
#num4=4
#
#u = TrialFunction(V)
#v = TestFunction(V)
#
#a=tr(dot((grad(u)).T, grad(v)))*dx
#L=-dot(normale,v)*ds(num1)-dot(normale,v)*ds(num2)-dot(normale,v)*ds(num3)-dot(normale,v)*ds(num4)
#
#u=Function(V)
#
#solve(a==L,u)
#
#H=assemble(grad(u)[0,0]*dx)
#M=assemble(grad(u)[1,1]*dx)
#A=assemble(grad(u)[0,1]*dx)
#C=assemble(grad(u)[1,0]*dx)
#print H,M,A,C
#ep=(1-pi*(0.25*r1**2+0.25*r2**2+0.25*r3**2+0.25*r4**2))
#D=((H/ep)+1)*ep
#print ("hi",r,ep,D)
#plot(u)
#plt.show()
#
###############################################
###################calcul erreur ##############
###############################################
#
#ured=Function(V)
#repertoire_read="Solution_homog_red_maill_init_multi/"
#lecture_solution_restart_hdf5(ured,repertoire_read,1,mesh)
#plot(ured)
#plt.show()
#
#num=Function(V)
#den=Function(V)
#
#u1=u.vector().array()
#u2=ured.vector().array()
#u3=u1-u2
#num.vector().set_local(u3)
#den.vector().set_local(u1)
#
#erreur=(sqrt(assemble(dot(num,num)*dx)))/(sqrt(assemble(dot(den,den)*dx)))
#
#print erreur











##############################################
########4 rectangle et cercle central#########
##############################################
##################test réduit#################
##############################################
#r1=0.05
#r2=0.1
#r3=0.15
#r4=0.2
##r1=r2=r3=r4

mesh = Mesh("Homog_base_maill_init_multi/mesh_complexe.xml.gz")



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
ep=(1-pi*(0.45**2)-4*0.12*0.12)
D=((H/ep)+1)*ep
print ("hi",ep,D)
plot(u)
plt.show()

##############################################
##################calcul erreur ##############
##############################################

ured=Function(V)
repertoire_read="Solution_homog_red_maill_init_multi/"
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











###############################################
#########4 rectangle et rectangle central######
###############################################
###################test réduit#################
###############################################
##r1=0.05
##r2=0.1
##r3=0.15
##r4=0.2
###r1=r2=r3=r4
#
#mesh = Mesh("Homog_base_maill_init_multi/mesh_complexe.xml.gz")
#
#
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
#ncel=V.dim()
#print ncel
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
#ds = Measure("ds")[boundaries]
#normale = FacetNormal(mesh) ## normale associe a chaque cellule
#
#num0=0
#num1=1
#num2=2
#num3=3
#num4=4
#
#u = TrialFunction(V)
#v = TestFunction(V)
#
#a=tr(dot((grad(u)).T, grad(v)))*dx
#L=-dot(normale,v)*ds(num0)-dot(normale,v)*ds(num1)-dot(normale,v)*ds(num2)-dot(normale,v)*ds(num3)-dot(normale,v)*ds(num4)
#
#u=Function(V)
#
#solve(a==L,u)
#
#H=assemble(grad(u)[0,0]*dx)
#M=assemble(grad(u)[1,1]*dx)
#A=assemble(grad(u)[0,1]*dx)
#C=assemble(grad(u)[1,0]*dx)
#print H,M,A,C
#ep=(1-0.3*0.3-4*0.15*0.15)
#D=((H/ep)+1)*ep
#print ("hi",ep,D)
#plot(u)
#plt.show()
#
###############################################
###################calcul erreur ##############
###############################################
#
#ured=Function(V)
#repertoire_read="Solution_homog_red_maill_init_multi/"
#lecture_solution_restart_hdf5(ured,repertoire_read,1,mesh)
#plot(ured)
#plt.show()
#
#num=Function(V)
#den=Function(V)
#
#u1=u.vector().array()
#u2=ured.vector().array()
#u3=u1-u2
#num.vector().set_local(u3)
#den.vector().set_local(u1)
#
#erreur=(sqrt(assemble(dot(num,num)*dx)))/(sqrt(assemble(dot(den,den)*dx)))
#
#print erreur
