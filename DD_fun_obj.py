from fenics import *

from dolfin import *
from mshr import *
import matplotlib.pyplot as plt
#import numpy as np
from math import sqrt
import sys



if __name__=='__main__':
 #
 tol=1e-10
 xinf=0.0
 yinf=0.0
 xsup=1.0
 ysup=1.0
 c_x=0.5
 c_y=0.5

### Une fonction, pour mailler un domaine avec une inclusion circulaire. Il en existe d'autres, voir homog_pod_multi ###

def raffinemment_maillage(cen,r,mesh):
 markers = CellFunction("bool", mesh)
 markers.set_all(False)
 c_x=cen[0]
 c_y=cen[1]
 for c in cells(mesh):
  # Mark cells with facet midpoints near y == 1.0
  for f in facets(c):
   if (sqrt((f.midpoint()[0]-c_x)**2+(f.midpoint()[1]-c_y)**2)<=1.2*r):
    markers[c] = True
  #new_mesh=refine(mesh, markers, redistribute=False)  
  new_mesh=refine(mesh, markers, redistribute=True)    
 return new_mesh

### Classes pour caractériser les bordures Gamma_sf et Gamma_ff , cas circulaire ###

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

# A renommer : interface #

class Obstacle(SubDomain):
 def inside(self, x, on_boundary):
  return (on_boundary and between((x[0]-c_x), (-r-tol, r+tol)) and between((x[1]-c_x), (-r-tol, r+tol)))

## Cas géométrique plus général, à réécrire ##

# Domaine d'intégration de chi, éventuellement pour des géométries plus compliquées que circulaire #

class Omega_fluide(SubDomain):
 def inside(self, x, on_boundary):
  return (((sqrt((x[0] - c_x)**2 + (x[1] - c_y)**2)-r) >= 0.0))








############################# Pour créer des maillages, avec des familles de cellules élémentaires #############################





def raffinement_maillage_cellule_centree(r,mesh):# Cellule centrée : un seul paramètre géométrique, le rayon de l'inclusion. 1.2*r<0.5 exigé.
 markers = CellFunction("bool", mesh)
 markers.set_all(False)
 c_x=0.5
 c_y=0.5
 for c in cells(mesh):
 # Mark cells with facet midpoints near y == 1.0
  for f in facets(c):
   if (sqrt((f.midpoint()[0]-c_x)**2+(f.midpoint()[1]-c_y)**2)<=1.2*r):
    markers[c] = True 
 mesh=refine(mesh, markers, redistribute=True) 
 return mesh

def raffinement_maillage_circ_per(cen,r,mesh):# Objectif : montrer que l'emplacement de l'inclusion périodique dans la cellule élémentaire ne change pas le coefficient de diffusion homogénéisé, calculé avec le tenseur khi
 markers = CellFunction("bool", mesh)
 markers.set_all(False)
 # on crée une liste des centres des inclusions voisines de la cellule élémentaire
 l_cen=[]
 for i in range(-1,2):
  for j in range(-1,2):
   l_cen.append([cen[0]+i,cen[1]+j])
 for c in cells(mesh):
 # Mark cells with facet midpoints near ...
  for f in facets(c):
   for cen_per in l_cen:
    if (sqrt((f.midpoint()[0]-cen_per[0])**2+(f.midpoint()[1]-cen_per[1])**2)<=1.2*r):
     markers[c] = True
 mesh=refine(mesh, markers, redistribute=True)
 return mesh

def creer_maill_circ(cen,r,res):#valable quel que soit la position de l'inclusion : centre, choisi aléatoirement. 1.2*r<0.5.
 if cen[0]+r*1.2<1 and cen[1]+r*1.2<1 and cen[0]-0.5*1.2>0 and cen[1]-0.5*1.2>0:#si l'inclusion est comprise dans la cellule
  rect=Rectangle(Point(0,0),Point(1,1))
  circle=Circle(Point(cen[0],cen[1]),r)
  domain=rect-circle
  mesh=generate_mesh(domain,res)
  # On raffine le long du bord de l'inclusion
  #mesh_aux=raffinement_maillage_cellule_centree(r,mesh)
  mesh_aux=raffinement_maillage_circ_per(cen,r,mesh)
  mesh=mesh_aux
 else:
  # Création de la cellule élémentaire avec inclusion
  rect=Rectangle(Point(-1,-1),Point(2,2))
  domain=rect
  #circle=Circle(Point(cen[0],cen[1]),r)
  l_cer=[Circle(Point(cen[0],cen[1]),r)]
  for i in range(-1,2):
   for j in range(-1,2):
    l_cer.append(Circle(Point(cen[0]+i,cen[1]+j),r))
  #domain_aux=rect-circle
  for cer_per in l_cer:
   domain=domain-cer_per
  domain=domain*Rectangle(Point(0,0),Point(1,1))
  # Création du permier maillage
  mesh=generate_mesh(domain,res)
  # On raffine le long du bord de l'inclusion
  mesh_aux=raffinement_maillage_circ_per(cen,r,mesh)
  mesh=mesh_aux
  #
 return(mesh)









