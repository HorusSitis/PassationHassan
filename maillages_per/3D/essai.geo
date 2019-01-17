
// Caracteristique maillage sphère
xc=0.5;
yc=0.5;
zc=0.5;
rayon=0.45;
pas_sphere=0.1;





/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
// Début du maillage de la sphère 
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
// Numero (moins 1) du premier point définissant la sphere
nps=700;
// Numero (moins 1) du premier cercle définissant la sphère
ncs=750;
// Numero (moins 1) de la première ligne définissant la sphère
nls=800;
// Numero (moins 1) de la première Surface définissant la sphère
nss=850;
// Numero de la surface totale définissant la sphere
num_surf_loop_sphere=900;
// Définition des points de construction de la sphère
Point(nps+1) = {xc,yc,zc,pas_sphere};
Point(nps+2) = {xc+rayon,yc,zc,pas_sphere};
Point(nps+3) = {xc,yc+rayon,zc,pas_sphere};
Point(nps+4) = {xc,yc,zc+rayon,pas_sphere};
Point(nps+5) = {xc-rayon,yc,zc,pas_sphere};
Point(nps+6) = {xc,yc-rayon,zc,pas_sphere};
Point(nps+7) = {xc,yc,zc-rayon,pas_sphere};

// Creation de 12 arcs de cercle
//Les arcs sur Gmsh sont définis à l'aide de la commande Circle. On utilise cette commande comme suit :
//Circle(numéro attribué à l'arc circulaire) = {numéro du point de départ, centre du cercle, point d'arrivée};
//On veut découper la sphère en 8 surfaces, pour ensuite considérer le volume enfermé à l'intérieur du regroupement ces surfaces.
Circle(nps+1) = {nps+2,nps+1,nps+3};
Circle(nps+2) = {nps+3,nps+1,nps+5};
