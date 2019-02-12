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

Point(npcy+1) = {xc,ymin,zc,pas_cube};
Point(npcy+2) = {xc+rayon,ymin,zc,pas_cube};
Point(npcy+3) = {xc-rayon,ymin,zc,pas_cube};

Point(npcy+4) = {xc,ymax,zc,pas_cube};
Point(npcy+5) = {xc+rayon,ymax,zc,pas_cube};
Point(npcy+6) = {xc-rayon,ymax,zc,pas_cube};

// Creation des contours du cylindre

Circle(nlcy+1)={npcy+2,npcy+1,npcy+3}
Circle(nlcy+2)={npcy+3,npcy+1,npcy+2}
Line(nlcy+3)={npcy+2,npcy+5}

Circle(nlcy+4)={npcy+5,npcy+4,npcy+6}
Circle(nlcy+5)={npcy+6,npcy+4,npcy+5}

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
Line Loop(nlcu+2) = {-(nlcu+11),nlcu+3,nlcu+9,nlcu+5};//,-(nlcy+1),-(nlcy+2)};
Line Loop(nlcy+6)={nlcy+1,nlcy+2}
Plane Surface(nPSc+2) = {nlcu+2,nlcy+7};
// Definition du contour fermé et de la surface située en ymax
Line Loop(nlcu+5) = {nlcu+12,nlcu+7,-(nlcu+10),nlcu+1};//,nlcy+4,nlcy+5};
Line Loop(nlcy+7)={nlcy+4,nlcy+5}
Plane Surface(nPSc+5) = {nlcu+5,nlcy+7};

/// Definition de la paroi du cylindre
Surface(num_surf_cylindre) = Extrude { { xc, 0, yc }, { 0, 1, 0 }, 2*Pi } {nlcy+3}

/// Périodicité : les orientations des lignes doivent correspondre

// // On impose la periodicité entre les surfaces d'équations xmin et xmax
Periodic Surface nPSc+4 {nlcu+11,-(nlcu+6), nlcu+12, -(nlcu+2)} = nPSc+3 {nlcu+9, nlcu+8, nlcu+10, nlcu+4};
// // On impose la periodicité entre les surfaces d'équations zmin et zmax
Periodic Surface nPSc+6 {-(nlcu+7),-(nlcu+6), -(nlcu+5), -(nlcu+8)} = nPSc+1 {nlcu+1, nlcu+2, nlcu+3, nlcu+4};

// On impose la periodicité entre les surfaces d'équations ymin et ymax
Periodic Surface nPSc+2 {nlcu+11,-(nlcu+3), nlcu+9, nlcu+5, nlcy+1, nlcy+2} = nPSc+5 {nlcu+9, nlcu+5, nlcu+11, nlcu+2, nlcy+4, nlcy+5};

// Surface du domaine fluide
Surface Loop(num_surf_loop_cube_trou) = {nPSc+1,-(nPSc+2),-(nPSc+3),-(nPSc+4),nPSc+5,nPSc+6};

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
// Creation du volume à mailler
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
Volume(num_vol) = {num_surf_loop_cube_trou,num_surf_cylindre};

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
