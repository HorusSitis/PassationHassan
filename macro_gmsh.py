### Pour utiliser la macro sur le terminal

import os
import sys

### Performances

import time

### Options

dimension=3
fig_todo='aff'

# Génération de maillages : apprentissage, fixe et test
appr=False
fixe=True
test=False

Nsnap=1

#res=1# résolution implicite, contenue dans les fichiers .geo
res=10
#res=20
#res=50
#res=100

## configurations en dimension 2
config='cer_un_som'#cer_un




config='sph_un'
if config=='sph_un':
 mesh_prefix='cubesphere_periodique_triangle'
##elif

if dimension==2:
 if config=='cer_un':
  mesh_prefix="maillage_trou2D"#"maillage_trou2d"#
 else:
  mesh_prefix="maillage_trou2D_som"
elif dimension==3:
 if config=='sph_un':
  mesh_prefix="cubesphere_periodique_triangle"
 elif config=='sph_un_som':
  mesh_prefix="cubesphere_periodique_triangle_som"
 #elif config=='cyl_ax':

dom_fixe="am"#""#"_0000"#"_0001"#

### On choisit un répertoire

#os.system("cd maillages_per/"+str(dimension)+"D")
os.chdir(os.getcwd() + "/maillages_per/"+str(dimension)+"D")

### Instructions à répéter en boucle : paramètre géométrique ou résolution variable, pour une configuration donnée
#geo_p=
#res_gmsh=

if appr:
 for n in range(1,1+Nsnap):
  r=n*0.05
  mesh_name=mesh_prefix+"_"+str(int(round(100*r,2)))
  if res!=100 and dimension==2:
   mesh_name=mesh_name+"_res"+str(res)
  ## res=100 : pas de suffixe
  if res!=10 and dimension==3:
   mesh_name=mesh_name+"_res"+str(res)
  ## res=10 : pas de suffixe
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
  #mesh_name=mesh_prefix+dom_fixe#+"_"
  mesh_name=mesh_prefix+"_am"+"_sur"+str(res)+"_fixe"
 if dimension==2:
  mesh_name="maillage_fixe2D_am"
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
 for r in [0.22,0.33,0.44]:
  mesh_name=mesh_prefix+"_"+str(int(round(100*r,2)))
  if res!=100 and dimension==2:
   mesh_name=mesh_name+"_res"+str(res)
  ## res=100 : pas de suffixe
  if res!=10 and dimension==3:
   mesh_name=mesh_name+"_res"+str(res)
  ## res=10 : pas de suffixe
  ###
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
