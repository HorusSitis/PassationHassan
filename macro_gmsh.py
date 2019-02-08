import os
import sys

### Options

dimension=3
fig_todo='aff'

#config='sph_cen'

#if config=='sph_cen':
# mesh_prefix=''
##elif

if dimension==2:
 mesh_prefix="maillage_trou2D"
elif dimension==3:
 mesh_prefix="cubesphere_periodique_triangle"

### On choisit un répertoire

#os.system("cd maillages_per/"+str(dimension)+"D")
os.chdir(os.getcwd() + "/maillages_per/"+str(dimension)+"D")

### Instructions à répéter en boucle : paramètre géométrique ou résolution variable, pour une configuration donnée
#geo_p=
#res_gmsh=

r=0.11# geo_p
mesh_name=mesh_prefix+"_"+str(int(round(100*r,2)))#+".xml"#"sur"+str(res_gmsh)+
print(mesh_name)

# Cas du domaine sans inclusion :
# mesh_name=mesh_prefix

## Génération d'un fichier .geo ? On commence avec un fichier unique et on modifie geo_p dans le code avant de sauvegarder sous le nom courant.




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
