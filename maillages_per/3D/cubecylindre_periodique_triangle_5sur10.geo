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
// Numero (moins 1) du premier point définissant le cylindre
npcy=200;

// Numero (moins 1) de la première ligne définissant le cube
nlcu=300;
// Numero (moins 1) de la première ligne définissant le cylindre
nlcy=400;

// Numero (moins 1) de la première Surface plane définissant le cube
nPSc=500;

// Nméro moins 1 de la première surface définissant le cylindre : quart de paroi
nCSc=550;

// Numero de la surface totale définissant le cube troué
num_surf_loop_cube_trou=600;
// Numero de la surface du cylindre
num_surf_cylindre=700;

// Numero du volume
num_vol=800;

// Numero des surfaces physiques
num_phys=900;

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

Point(npcy+1) = {xc,ymin,zc,pas_cylindre};
Point(npcy+2) = {xc+rayon,ymin,zc,pas_cylindre};
Point(npcy+3) = {xc,ymin,zc+rayon,pas_cylindre};
Point(npcy+4) = {xc-rayon,ymin,zc,pas_cylindre};
Point(npcy+5) = {xc,ymin,zc-rayon,pas_cylindre};

Point(npcy+6) = {xc,ymax,zc,pas_cylindre};
Point(npcy+7) = {xc+rayon,ymax,zc,pas_cylindre};
Point(npcy+8) = {xc,ymax,zc+rayon,pas_cylindre};
Point(npcy+9) = {xc-rayon,ymax,zc,pas_cylindre};
Point(npcy+10) = {xc,ymax,zc-rayon,pas_cylindre};


// Creation des contours du cylindre

Circle(nlcy+1)={npcy+2,npcy+1,npcy+3};
Circle(nlcy+2)={npcy+3,npcy+1,npcy+4};
Circle(nlcy+3)={npcy+4,npcy+1,npcy+5};
Circle(nlcy+4)={npcy+5,npcy+1,npcy+2};

Line(nlcy+5)={npcy+2,npcy+7};

Circle(nlcy+6)={npcy+7,npcy+6,npcy+8};
Circle(nlcy+7)={npcy+8,npcy+6,npcy+9};
Circle(nlcy+8)={npcy+9,npcy+6,npcy+10};
Circle(nlcy+9)={npcy+10,npcy+6,npcy+7};

/// Ajouts pour la surface du cylindre
Line(nlcy+12)={npcy+3,npcy+8};
Line(nlcy+13)={npcy+4,npcy+9};
Line(nlcy+14)={npcy+5,npcy+10};
Line Loop(nlcy+15)={nlcy+1,nlcy+12,-(nlcy+6),-(nlcy+5)};
Line Loop(nlcy+16)={nlcy+2,nlcy+13,-(nlcy+7),-(nlcy+12)};
Line Loop(nlcy+17)={nlcy+3,nlcy+14,-(nlcy+8),-(nlcy+13)};
Line Loop(nlcy+18)={nlcy+4,nlcy+5,-(nlcy+9),-(nlcy+14)};

//Curve(nl)


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

/// Definition des contours des faces trouées : on ajoute les parois des trous
// Definition du contour fermé et de la surface située en ymin
Line Loop(nlcu+2) = {-(nlcu+11),nlcu+3,nlcu+9,nlcu+5};
Line Loop(nlcy+10)={nlcy+1,nlcy+2,nlcy+3,nlcy+4};
Plane Surface(nPSc+2) = {nlcu+2,nlcy+10};
// Definition du contour fermé et de la surface située en ymax
Line Loop(nlcu+5) = {nlcu+12,nlcu+7,-(nlcu+10),nlcu+1};
Line Loop(nlcy+11)={nlcy+6,nlcy+7,nlcy+8,nlcy+9};
Plane Surface(nPSc+5) = {nlcu+5,nlcy+11};

/// Definition de la paroi du cylindre
//num_surf_cylindre = 
//Extrude {{0,1,0},{xc,0,zc},2*Pi} {Curve{nlcy+5};};
//out[]
//rev[] = 

Ruled Surface(nCSc+1)={nlcy+15};
Ruled Surface(nCSc+2)={nlcy+16};
Ruled Surface(nCSc+3)={nlcy+17};
Ruled Surface(nCSc+4)={nlcy+18};
Surface Loop(num_surf_cylindre) = {nCSc+1,nCSc+2,nCSc+3,nCSc+4};

////num_surf_cylindre = Extrude{{0,1,0},{xc,0,zc},2*Pi}{ Curve{nlcy+5}; };
//Ruled 
//Surface Loop(700) = {rev[]};
//Surface num_surf_cylindre {rev[1]};
//Surface(num_surf_cylindre) = Surface{rev[1]};

/// Périodicité : les orientations des lignes doivent correspondre

// // On impose la periodicité entre les surfaces d'équations xmin et xmax
Periodic Surface nPSc+3 {nlcu+11,-(nlcu+6), nlcu+12, -(nlcu+2)} = nPSc+4 {nlcu+9, nlcu+8, nlcu+10, nlcu+4};
// // On impose la periodicité entre les surfaces d'équations zmin et zmax
Periodic Surface nPSc+1 {-(nlcu+7),-(nlcu+6), -(nlcu+5), -(nlcu+8)} = nPSc+6 {nlcu+1, nlcu+2, nlcu+3, nlcu+4};

// On impose la periodicité entre les surfaces d'équations ymin et ymax
Periodic Surface nPSc+2 {nlcu+11,-(nlcu+5), -(nlcu+9), -(nlcu+3), nlcy+1, nlcy+2, nlcy+3, nlcy+4} = nPSc+5 {nlcu+12, nlcu+7, -(nlcu+10), nlcu+1, nlcy+6, nlcy+7, nlcy+8, nlcy+9};

// Surface du domaine fluide
Surface Loop(num_surf_loop_cube_trou) = {nPSc+1,-(nPSc+2),-(nPSc+3),-(nPSc+4),nPSc+5,nPSc+6};

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
// Creation du volume à mailler
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
Volume(num_vol) = {num_surf_loop_cube_trou,num_surf_cylindre};//Surface{rev[1]}};//

//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
// Definition des surfaces physiques : numérotation pour les conditions aux limites //
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////

Physical Surface(num_phys+1) = {nPSc+6}; //surface située en zmin
Physical Surface(num_phys+2) = {nPSc+1}; //surface située en zmax
Physical Surface(num_phys+3) = {nPSc+4}; // surface située en xmin
Physical Surface(num_phys+4) = {nPSc+3}; // surface située en ymax
Physical Surface(num_phys+5) = {nPSc+2}; // surface située en xmax
Physical Surface(num_phys+6) = {nPSc+5}; // surface située en ymin

Physical Surface(num_phys+7) = {num_surf_cylindre}; // surface de la sphère

/////////////////////////////////////////////////////////////////
// Definition du volume physique
/////////////////////////////////////////////////////////////////
Physical Volume(num_phys+8) = {num_vol};
