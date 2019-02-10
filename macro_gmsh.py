import os
import sys

### Options

dimension=2
fig_todo='aff'

# Génération de maillages : apprentissage, fixe et test
appr=False
fixe=False
test=True

Nsnap=1

config='cer_un_som'#cer_un

res=1# valeur par défaut, jamais utilisée
#res=20

#config='sph_cen'
#if config=='sph_cen':
# mesh_prefix=''
##elif

if dimension==2:
 if config=='cer_un':
  mesh_prefix="maillage_trou2D"#"maillage_trou2d"#
 else:
  mesh_prefix="maillage_trou2D_som"
elif dimension==3:
 mesh_prefix="cubesphere_periodique_triangle"

dom_fixe="_0001"#""#"_0000"#"_0001"#

### On choisit un répertoire

#os.system("cd maillages_per/"+str(dimension)+"D")
os.chdir(os.getcwd() + "/maillages_per/"+str(dimension)+"D")

### Instructions à répéter en boucle : paramètre géométrique ou résolution variable, pour une configuration donnée
#geo_p=
#res_gmsh=

if appr:
 for n in range(1,1+Nsnap):
  r=n*0.05
  #if r==0.05:
  # mesh_name=mesh_prefix+"_0"+str(int(round(100*r,2)))
  #else:
  mesh_name=mesh_prefix+"_"+str(int(round(100*r,2)))
  if res==20:
   mesh_name=mesh_name+"_res"+str(res)
  ## res=100 : pas de suffixe
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
  mesh_name=mesh_prefix+dom_fixe#+"_"
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
  if res==20:
   mesh_name=mesh_name+"_res"+str(res)
  ## res=100 : pas de suffixe
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

sys.exit("Création de maillages périodiques terminée")
