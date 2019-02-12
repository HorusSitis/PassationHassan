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
npcu=100;
// Numero (moins 1) du premier point définissant les sections frontales du cylindre
npcy=200;

// Numero (moins 1) de la première ligne définissant le cube
nlcu=300;
// Numero (moins 1) de la première ligne définissant le cylindre
nlcy=400;

// Numero (moins 1) de la première Surface plane définissant le cube
nPSc=500;
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
Point(npcu+1) = {xmin,ymin,zmin,pas_cube};
Point(npcu+2) = {xmax,ymin,zmin,pas_cube};
Point(npcu+3) = {xmax,ymax,zmin,pas_cube};
Point(npcu+4) = {xmin,ymax,zmin,pas_cube};
Point(npcu+5) = {xmin,ymin,zmax,pas_cube};
Point(npcu+6) = {xmax,ymin,zmax,pas_cube};
Point(npcu+7) = {xmax,ymax,zmax,pas_cube};
Point(npcu+8) = {xmin,ymax,zmax,pas_cube};

// Creation des 12 arêtes du cube
Line(nlcu+1) = {npcu+4,npcu+3};
Line(nlcu+2) = {npcu+3,npcu+2};
Line(nlcu+3) = {npcu+2,npcu+1};
Line(nlcu+4) = {npcu+1,npcu+4};

Line(nlcu+5) = {npcu+5,npcu+6};
Line(nlcu+6) = {npcu+6,npcu+7};
Line(nlcu+7) = {npcu+7,npcu+8};
Line(nlcu+8) = {npcu+8,npcu+5};

Line(nlcu+9) = {npcu+1,npcu+5};
Line(nlcu+10) = {npcu+4,npcu+8};
Line(nlcu+11) = {npcu+2,npcu+6};
Line(nlcu+12) = {npcu+3,npcu+7};

// Creation des points d'intersection entre l'axe du cylindre et les faces frontales

Point(npcy+1) = {xc,ymin,zc,pas_cube};
Point(npcy+2) = {xc+rayon,ymin,zc,pas_cube};
Point(npcy+3) = {xc-rayon,ymin,zc,pas_cube};

Point(npcy+4) = {xc,ymax,zc,pas_cube};
Point(npcy+5) = {xc+rayon,ymax,zc,pas_cube};
Point(npcy+6) = {xc-rayon,ymax,zc,pas_cube};

// Creation des sections du cylindre sur les faces frontales

Circle(nlcy+1)={npcy+2,npcy+1,npcy+3}
Circle(nlcy+2)={npcy+5,npcy+4,npcy+6}

/// Definition des contours carrés : numéros de dé

// Definition du contour fermé et de la surface située en zmin
Line Loop(nlcu+6) = {nlcu+3,nlcu+4,nlcu+1,nlcu+2};
Plane Surface(nPSc+6) = {nlcu+6};
// Definition du contour fermé et de la surface située en zmax
Line Loop(nlcu+1) = {nlcu+5,nlcu+6,nlcu+7,nlcu+8};
Plane Surface(nPSc+1) = {nlcu+1};
// Definition du contour fermé et de la surface située en xmin
Line Loop(nlcu+4) = {nlcu+9,-(nlcu+8),-(nlcu+10),-(nlcu+4)};
Plane Surface(nPSc+4) = {nlcu+4};
// Definition du contour fermé et de la surface située en xmax
Line Loop(nlcu+3) = {nlcu+11,nlcu+6,-(nlcu+12),nlcu+2};
Plane Surface(nPSc+3) = {nlcu+3};

/// Definition des contours des faces trouées




// Definition du contour fermé et de la surface située en ymin
Line Loop(nlc+2) = {-(nlcu+11),nlcu+5,nlcu+9,nlcu+3,-(nlcy+1)};
Plane Surface(nPSc+2) = {nlc+2};
// Definition du contour fermé et de la surface située en ymax
Line Loop(nlc+5) = {nlc+12,nlc+7,-(nlc+10),nlc+1,nlcy+2};
Plane Surface(nPSc+5) = {nlc+5};

/// Definition des parois du cylindre



/// Périodicité : les orientations des lignes doivent correspondre

// // On impose la periodicité entre les surfaces d'équations xmin et xmax
Periodic Surface nPSc+4 {nlc+11,-(nlc+6), nlc+12, -(nlc+2)} = nPSc+3 {nlc+9, nlc+8, nlc+10, nlc+4};
// // On impose la periodicité entre les surfaces d'équations zmin et zmax
Periodic Surface nPSc+6 {-(nlc+7),-(nlc+6), -(nlc+5), -(nlc+8)} = nPSc+1 {nlc+1, nlc+2, nlc+3, nlc+4};

// On impose la periodicité entre les surfaces d'équations ymin et ymax
Periodic Surface nPSc+2 {nlc+10,-(nlc+7), nlc+12, -(nlc+1)} = nPSc+5 {nlc+9, nlc+5, nlc+11, nlc+2};

// Surface du domaine fluide
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
// Creation du volume à mailler
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
Volume(num_vol) = {num_surf_loop_cube,num_surf_loop_sphere};


//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
// Definition des surfaces physiques : numérotation pour les conditions aux limites //
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////

Physical Surface(num_phys+1) = {nPSc+15}; //surface située en zmin
Physical Surface(num_phys+2) = {nPSc+17}; //surface située en zmax
Physical Surface(num_phys+3) = {nPSc+19}; // surface située en xmin
Physical Surface(num_phys+4) = {nPSc+21}; // surface située en ymax
Physical Surface(num_phys+5) = {nPSc+23}; // surface située en xmax
Physical Surface(num_phys+6) = {nPSc+25}; // surface située en ymin

Physical Surface(num_phys+7) = {num_surf_loop_sphere}; // surface de la sphère

/////////////////////////////////////////////////////////////////
// Definition du volume physique
/////////////////////////////////////////////////////////////////
Physical Volume(num_phys+3) = {num_vol};
