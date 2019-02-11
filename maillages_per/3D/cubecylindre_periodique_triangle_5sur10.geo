// Caracteristique maillage cube 
pas_cube=0.1;// pas_cube : pas pour mailler le cube

xmin=0.;
xmax=1.;
ymin=0.;
ymax=1.;
zmin=0.;
zmax=1.;

// Caracteristique maillage cylindre
xc=0.5;

zc=0.5;

rayon=0.05;
pas_cylindre=0.1;


///

// Numero (moins 1) du premier point définissant le cube
npc=100;
// Numero (moins 1) de la première ligne définissant le cube
nlc=200;
// Numero (moins 1) de la première Surface plane définissant le cube
nPSc=300;
// Numero (moins 1) de la première Surface définissant le cube
nsc=400;
// Numero de la surface totale définissant le cube
num_surf_loop_cube=600;
// Numero du volume
num_vol=500;
// Numero des surfaces physiques
num_phys=10;

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
// Début du maillage 
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

// Creation des 8 points du bord du cube
Point(npc+1) = {xmin,ymin,zmin,pas_cube};
Point(npc+2) = {xmax,ymin,zmin,pas_cube};
Point(npc+3) = {xmax,ymax,zmin,pas_cube};
Point(npc+4) = {xmin,ymax,zmin,pas_cube};
Point(npc+5) = {xmin,ymin,zmax,pas_cube};
Point(npc+6) = {xmax,ymin,zmax,pas_cube};
Point(npc+7) = {xmax,ymax,zmax,pas_cube};
Point(npc+8) = {xmin,ymax,zmax,pas_cube};

// Creation des 8 arêtes du cube
Line(nlc+1) = {npc+4,npc+3};
Line(nlc+2) = {npc+3,npc+2};
Line(nlc+3) = {npc+2,npc+1};
Line(nlc+4) = {npc+1,npc+4};
Line(nlc+6) = {npc+5,npc+6};
Line(nlc+7) = {npc+6,npc+7};
Line(nlc+8) = {npc+7,npc+8};
Line(nlc+9) = {npc+8,npc+5};
Line(nlc+10) = {npc+1,npc+5};
Line(nlc+11) = {npc+4,npc+8};
Line(nlc+12) = {npc+2,npc+6};
Line(nlc+13) = {npc+3,npc+7};

// Creation des points d'inersection entre l'axe du cylindre et les faces frontales


// Creation des sections du cylindre sur les faces frontales



// Definition du contour fermé et de la surface située en zmin
Line Loop(nlc+14) = {nlc+3,nlc+4,nlc+1,nlc+2};
Plane Surface(nPSc+15) = {nlc+14};
// Definition du contour fermé et de la surface située en zmax
Line Loop(nlc+16) = {nlc+6,nlc+7,nlc+8,nlc+9};
Plane Surface(nPSc+17) = {nlc+16};
// Definition du contour fermé et de la surface située en xmin
Line Loop(nlc+18) = {nlc+10,-(nlc+9),-(nlc+11),-(nlc+4)};
Plane Surface(nPSc+19) = {nlc+18};
// Definition du contour fermé et de la surface située en xmax
Line Loop(nlc+22) = {nlc+12,nlc+7,-(nlc+13),nlc+2};
Plane Surface(nPSc+23) = {nlc+22};



// Definition des parois du cylindre


// // On impose la periodicité entre les surfaces d'équations xmin et xmax
Periodic Surface nPSc+23 {nlc+12,-(nlc+7), nlc+13, -(nlc+2)} = nPSc+19 {nlc+10, nlc+9, nlc+11, nlc+4};
// On impose la periodicité entre les surfaces d'équations ymin et ymax
Periodic Surface nPSc+21 {nlc+11,-(nlc+8), nlc+13, -(nlc+1)} = nPSc+25 {nlc+10, nlc+6, nlc+12, nlc+3};
// // On impose la periodicité entre les surfaces d'équations zmin et zmax
Periodic Surface nPSc+17 {-(nlc+8),-(nlc+7), -(nlc+6), -(nlc+9)} = nPSc+15 {nlc+1, nlc+2, nlc+3, nlc+4};


Surface Loop(num_surf_loop_cube) = {nPSc+17,-(nPSc+25),-(nPSc+23),-(nPSc+21),nPSc+19,nPSc+15};

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
// Fin du maillage du cube 
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////


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
Circle(nps+3) = {nps+5,nps+1,nps+6};
Circle(nps+4) = {nps+6,nps+1,nps+2};
Circle(nps+5) = {nps+2,nps+1,nps+7};
Circle(nps+6) = {nps+7,nps+1,nps+5};
Circle(nps+7) = {nps+5,nps+1,nps+4};
Circle(nps+8) = {nps+4,nps+1,nps+2};
Circle(nps+9) = {nps+6,nps+1,nps+7};
Circle(nps+10) = {nps+7,nps+1,nps+3};
Circle(nps+11) = {nps+3,nps+1,nps+4};
Circle(nps+12) ={nps+4,nps+1,nps+6};


//Définition des contour de la sphère avec la commande Line Loop.
// Line Loop (numéro du contour attribué au contours) = {arcs ou lignes à relier};
//Avec cette commande, il faut faire attention à l'orientation des lignes et des arcs. 
// Pour changer l'orientation d'une ligne ou d'un arc on rajoute le signe - avant son numéro.
Line Loop(nls+1) = {nps+1,nps+11,nps+8};
Line Loop(nls+2) = {nps+2,nps+7,-(nps+11)};
Line Loop(nls+3) = {nps+3,-(nps+12),-(nps+7)};
Line Loop(nls+4) = {nps+4,-(nps+8),nps+12};
Line Loop(nls+5) = {nps+5,nps+10,-(nps+1)};
Line Loop(nls+6) = {-(nps+2),-(nps+10),nps+6};
Line Loop(nls+7) = {-(nps+3),-(nps+6),-(nps+9)};
Line Loop(nls+8) = {-(nps+4),nps+9,-(nps+5)};


// 5/ Surfaces :
//Surface : permet de définir des surfaces sphériques à partir des contours.
Ruled Surface(nss+1) = {nls+1};
Ruled Surface(nss+2) = {nls+2};
Ruled Surface(nss+3) = {nls+3};
Ruled Surface(nss+4) = {nls+4};
Ruled Surface(nss+5) = {nls+5};
Ruled Surface(nss+6) = {nls+6};
Ruled Surface(nss+7) = {nls+7};
Ruled Surface(nss+8) = {nls+8};

//Et Surface Loop : permet de définir la surface fermée engendrée par les 8 surfaces définie ci-dessous.
Surface Loop (num_surf_loop_sphere) = {nss+1,nss+2,nss+3,nss+4,nss+5,nss+6,nss+7,nss+8};


/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
// Fin du maillage du cube 
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
// Creation du bolume à mailler
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
Volume(num_vol) = {num_surf_loop_cube,num_surf_loop_sphere};


/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
// Definition des surfaces physiques
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
Physical Surface(num_phys+1) = {nPSc+15}; //surface située en zmin
Physical Surface(num_phys+2) = {nPSc+17}; //surface située en zmax
Physical Surface(num_phys+3) = {nPSc+19}; // surface située en xmin
Physical Surface(num_phys+4) = {nPSc+21}; // surface située en ymax
Physical Surface(num_phys+5) = {nPSc+23}; // surface située en xmax
Physical Surface(num_phys+6) = {nPSc+25}; // surface située en ymin
Physical Surface(num_phys+7) = {num_surf_loop_sphere}; // surface de la sphère
Physical Volume(num_phys+3) = {num_vol}; 
