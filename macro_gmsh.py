### Pour utiliser la macro sur le terminal

import os
import sys

### Performances

import time

### Options

dimension=2
fig_todo=''

# Génération de maillages : apprentissage, fixe et test
appr=False
fixe=True
test=False


Nsnap=8

res_name=False
res=10
#res=20
#res=50
#res=100

## configurations en dimension 2
if dimension==2:
 #config='cer_un'
 #config='cer_un_som'
 config='compl'
 #geo_p='diag'
 geo_p='hor'

## configurations en dimension 3
if dimension==3:
 #config='sph_un'
 config='cyl_un'

if dimension==2:
 if config=='cer_un':
  mesh_prefix="maillage_trou2D"#"maillage_trou2d"#
 elif config=='cer_un_som':
  mesh_prefix="maillage_trou2D_som"
 elif config=='compl':
  mesh_prefix="maillage_trous2D_"+geo_p  
elif dimension==3:
 if config=='sph_un':
  mesh_prefix="cubesphere_periodique_triangle"
 elif config=='sph_un_som':
  mesh_prefix="cubesphere_periodique_triangle_som"
 elif config=='cyl_un':
  mesh_prefix='cubecylindre_periodique_triangle'

dom_fixe="am"

### On choisit un répertoire

#os.system("cd maillages_per/"+str(dimension)+"D")
os.chdir(os.getcwd() + "/maillages_per/"+str(dimension)+"D")

### Instructions à répéter en boucle : paramètre géométrique ou résolution variable, pour une configuration donnée
#geo_p=
#res_gmsh=

if appr:
 for n in range(1,2+Nsnap):
  if dimension==2 and geo_p=='hor':
   r=n*0.04+0.01
  elif dimension==3 and geo_p=='lat':
   r=n*0.04+0.01##??
  else:
   r=n*0.05
  mesh_name=mesh_prefix+"_"+str(int(round(100*r,2)))
  if res_name and dimension==2:
   mesh_name=mesh_name+"_res"+str(res)
  ## res=100 : pas de suffixe
  if res_name and dimension==3:
   mesh_name=mesh_name+"sur"+str(res)
  ## res=10 : pas de suffixe dans le cas sphérique, voir avec res_name
  ###
  ## Génération d'un fichier .geo ? On commence avec un fichier unique et on modifie geo_p dans le code avant de sauvegarder sous le nom courant.
  print(mesh_name)
  ## Visualisation du fichier .geo
  print("gmsh "+mesh_name+".geo")
  os.system("gmsh "+mesh_name+".geo")
  ## Conversion en .msh
  os.system("gmsh -"+str(dimension)+" "+mesh_name+".geo")
  ## Affichage du maillage obtenu
  if fig_todo=='aff':
   os.system("gmsh "+mesh_name+".msh")
  ## Conversion en .xml avec dolfin pour FEniCS
  os.system("dolfin-convert "+mesh_name+".msh "+mesh_name+".xml")
elif fixe:
 if dimension==3:
  if dom_fixe=="am":
   mesh_name=mesh_prefix+"_am"+"_sur"+str(res)+"_fixe"
  elif dom_fixe=="0001":
   mesh_name=mesh_prefix+"_sur"+str(res)+"_"+dom_fixe+"fixe"
 if dimension==2:
  if config!='compl':
   mesh_name="maillage_fixe2D_am"
  else:
   mesh_name=mesh_prefix+"_fixe"
   print(mesh_name)
 ## Génération d'un fichier .geo ? On commence avec un fichier unique et on modifie geo_p dans le code avant de sauvegarder sous le nom courant.
 print(mesh_name)
 ## Visualisation du fichier .geo
 print("gmsh "+mesh_name+".geo")
 os.system("gmsh "+mesh_name+".geo")
 ## Conversion en .msh
 os.system("gmsh -"+str(dimension)+" "+mesh_name+".geo")
 ## Affichage du maillage obtenu
 if fig_todo=='aff':
  os.system("gmsh "+mesh_name+".msh")
 ## Conversion en .xml avec dolfin pour FEniCS
 os.system("dolfin-convert "+mesh_name+".msh "+mesh_name+".xml")
elif test:
 for r in [0.11,0.22,0.33,0.44]:# pour compl
  mesh_name=mesh_prefix+"_"+str(int(round(100*r,2)))
  if res_name and dimension==2:
   mesh_name=mesh_name+"_res"+str(res)
  ## res=100 : pas de suffixe
  if res_name and dimension==3:
   mesh_name=mesh_name+"sur"+str(res)
  ## res=10 : pas de suffixe dans le cas sphérique, voir avec res_name
  else:
   mesh_name=mesh_prefix+"_"+str(int(round(100*r,2)))
  ## Génération d'un fichier .geo ? On commence avec un fichier unique et on modifie geo_p dans le code avant de sauvegarder sous le nom courant.
  print(mesh_name)
  ## Visualisation du fichier .geo
  print("gmsh "+mesh_name+".geo")
  os.system("gmsh "+mesh_name+".geo")
  ## Conversion en .msh
  start=time.time()
  os.system("gmsh -"+str(dimension)+" "+mesh_name+".geo")
  end=time.time()
  print("temps de génération du maillage : "+str(end-start)+" secondes")
  ## Affichage du maillage obtenu
  if fig_todo=='aff':
   os.system("gmsh "+mesh_name+".msh")
  ## Conversion en .xml avec dolfin pour FEniCS
  start=time.time()
  os.system("dolfin-convert "+mesh_name+".msh "+mesh_name+".xml")
  end=time.time()
  print("temps de conversion du maillage : "+str(end-start)+" secondes")

sys.exit("Création de maillages périodiques terminée")
