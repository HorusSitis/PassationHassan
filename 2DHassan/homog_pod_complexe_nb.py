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

R_dim=len(R)

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

#domaine_fixe=Rectangle(Point(xinf,yinf),Point(xsup,ysup))
#mesh_fixe=generate_mesh(domaine_fixe,100)
mesh_fixe = Mesh("Solutions_homog_interp_complexe_nb/mesh_complexe_nb.xml.gz")


V_fixe=VectorFunctionSpace(mesh_fixe, "P", 2, constrained_domain=PeriodicBoundary())

nb_noeuds = V_fixe.dim()
print nb_noeuds

Usnap=np.zeros((nb_noeuds,R_dim))

repertoire_read="Solutions_homog_interp_complexe_nb/"

u = Function(V_fixe)

for n in range(R_dim):
      r=R[n]
      lecture_fichiers_restart_hdf5(u,repertoire_read,r,mesh_fixe) 
      Usnap[:,n]=u.vector().array()

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

nb_modes=18

print nb_modes

#choix de la base POD
G=np.zeros((nb_noeuds,nb_modes))
#base=Function(V_fixe)
for i in range(nb_modes):
    G[:,i]=phi[:,i]
#    base.vector().set_local(G[:,i])
#    plot(base)
#    plt.show()
#    print base(0.5,0.5),base(0.1,0.99),base(0.2,0.9),base(0.4,0.03)

################################################################
########Interpolation de la base dans chaque maillage###########
################################################################


r=0.123

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

repertoire_final_1="Homog_base_maill_init_complexe_nb/"

File("Homog_base_maill_init_complexe_nb/mesh_complexe_nb.xml.gz") << mesh

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

repertoire_final_2="Solution_homog_red_maill_init_complexe_nb/"

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
ep=(1-pi*((r**2)+0.45**2))
D=((H/ep)+1)*ep
print ("hi",r,ep,D)



plot(urec_red)
plt.show()